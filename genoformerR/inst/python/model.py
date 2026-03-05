# genoformerR/inst/python/model.py  (v0.2.0)
# GenoFormer PyTorch model with dual-level ancestry conditioning:
#   - Global ancestry: FiLM at every transformer block (genome-wide background)
#   - Local ancestry:  Per-SNP token embedding (haplotype-mosaic context)
# Sourced by gf_init() via reticulate::py_run_file()

import torch
import torch.nn as nn
import torch.nn.functional as F


class AncestryFiLM(nn.Module):
    """FiLM conditioning on global ancestry — applied to each block residual."""
    def __init__(self, anc_dim: int, hidden_dim: int):
        super().__init__()
        self.gamma = nn.Linear(anc_dim, hidden_dim)
        self.beta  = nn.Linear(anc_dim, hidden_dim)
        nn.init.ones_(self.gamma.weight);  nn.init.zeros_(self.gamma.bias)
        nn.init.zeros_(self.beta.weight);  nn.init.zeros_(self.beta.bias)

    def forward(self, x, global_anc):
        g = self.gamma(global_anc).unsqueeze(1)
        b = self.beta(global_anc).unsqueeze(1)
        return g * x + b


class LocalAncestryProjector(nn.Module):
    """Per-SNP local ancestry projection: (B, L, K) -> (B, L, d_local)."""
    def __init__(self, k_pops: int, d_local: int):
        super().__init__()
        self.proj = nn.Sequential(
            nn.Linear(k_pops, d_local),
            nn.GELU(),
            nn.Linear(d_local, d_local)
        )

    def forward(self, local_anc):
        return self.proj(local_anc)


class SNPEmbedding(nn.Module):
    """
    Token embedding for each SNP position.
    v0.2 adds local ancestry fusion: token = f(dosage, log|b|, LD, chr, local_anc)
    """
    def __init__(self, n_ld_blocks, n_chrom, d_model, k_pops=0, d_local=0):
        super().__init__()
        quarter = d_model // 4
        self.dosage_proj = nn.Linear(1, quarter)
        self.beta_proj   = nn.Linear(1, quarter)
        self.ld_emb      = nn.Embedding(n_ld_blocks + 1, quarter, padding_idx=0)
        self.chr_emb     = nn.Embedding(n_chrom + 1,     quarter, padding_idx=0)

        self.use_local_anc = k_pops > 0 and d_local > 0
        if self.use_local_anc:
            self.local_proj = LocalAncestryProjector(k_pops, d_local)
            self.fusion = nn.Sequential(
                nn.Linear(d_model + d_local, d_model),
                nn.LayerNorm(d_model)
            )
        else:
            self.norm = nn.LayerNorm(d_model)

    def forward(self, dosage, log_beta, ld_block, chrom, local_anc=None):
        base = torch.cat([
            self.dosage_proj(dosage),
            self.beta_proj(log_beta),
            self.ld_emb(ld_block),
            self.chr_emb(chrom)
        ], dim=-1)
        if self.use_local_anc and local_anc is not None:
            la = self.local_proj(local_anc)
            return self.fusion(torch.cat([base, la], dim=-1))
        return self.norm(base)


class GenoFormerBlock(nn.Module):
    """Transformer block; global ancestry FiLM on attention residual."""
    def __init__(self, d_model, n_heads, anc_dim, dropout=0.1):
        super().__init__()
        self.attn  = nn.MultiheadAttention(d_model, n_heads,
                                            dropout=dropout, batch_first=True)
        self.film  = AncestryFiLM(anc_dim, d_model)
        self.ff    = nn.Sequential(
            nn.Linear(d_model, d_model * 4), nn.GELU(),
            nn.Dropout(dropout), nn.Linear(d_model * 4, d_model)
        )
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)

    def forward(self, x, global_anc, attn_mask=None):
        attn_out, _ = self.attn(x, x, x, attn_mask=attn_mask)
        x = self.norm1(x + self.film(attn_out, global_anc))
        x = self.norm2(x + self.ff(x))
        return x


class GenoFormer(nn.Module):
    """
    GenoFormer v0.2 — dual-level ancestry-aware PRS transformer.

    Dual conditioning:
      Local  ancestry: fused into each SNP token  (per-position, K-vector)
      Global ancestry: FiLM on every block output (genome-wide background)

    Falls back to global-only (v0.1 behaviour) when use_local_anc=False
    or local_anc=None is passed to forward().
    """
    def __init__(self,
                 n_snps=80_000, n_ld_blocks=1_700, n_chrom=22,
                 d_model=256, n_heads=8, n_layers=6,
                 anc_dim=5, k_pops=5, d_local=64,
                 dropout=0.1, use_local_anc=True):
        super().__init__()
        self.use_local_anc = use_local_anc
        self.embed = SNPEmbedding(
            n_ld_blocks, n_chrom, d_model,
            k_pops  = k_pops  if use_local_anc else 0,
            d_local = d_local if use_local_anc else 0
        )
        self.blocks    = nn.ModuleList([
            GenoFormerBlock(d_model, n_heads, anc_dim, dropout)
            for _ in range(n_layers)
        ])
        self.cls_token = nn.Parameter(torch.randn(1, 1, d_model) * 0.02)
        self.prs_head  = nn.Sequential(
            nn.Linear(d_model, 64), nn.GELU(),
            nn.Dropout(dropout), nn.Linear(64, 1)
        )
        self.pop_head  = nn.Linear(d_model, anc_dim)

    def forward(self, dosage, log_beta, ld_block, chrom,
                global_anc, local_anc=None):
        B = dosage.shape[0]
        x = self.embed(
            dosage.unsqueeze(-1), log_beta.unsqueeze(-1),
            ld_block, chrom,
            local_anc if self.use_local_anc else None
        )
        cls = self.cls_token.expand(B, -1, -1)
        x   = torch.cat([cls, x], dim=1)
        for block in self.blocks:
            x = block(x, global_anc)
        cls_repr = x[:, 0]
        return self.prs_head(cls_repr).squeeze(-1), self.pop_head(cls_repr)


class GenoFormerPhased(nn.Module):
    """
    Haplotype-aware variant — processes hap0 and hap1 separately through
    shared GenoFormer weights. PRS = hap0_prs + hap1_prs (additive diploid).
    Requires phased dosage (binary 0/1) and per-haplotype local ancestry.
    """
    def __init__(self, **kwargs):
        super().__init__()
        self.base = GenoFormer(**kwargs)

    def forward(self, hap0_dos, hap1_dos, log_beta, ld_block, chrom,
                global_anc, hap0_local_anc, hap1_local_anc):
        prs0, pop0 = self.base(hap0_dos, log_beta, ld_block, chrom,
                               global_anc, hap0_local_anc)
        prs1, pop1 = self.base(hap1_dos, log_beta, ld_block, chrom,
                               global_anc, hap1_local_anc)
        return prs0 + prs1, (pop0 + pop1) / 2.0
