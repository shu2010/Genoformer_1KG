# genoformerR — Function Reference (v0.3)

## Overview

genoformerR implements a transformer-based, ancestry-aware polygenic risk scoring
pipeline. Functions are grouped into seven stages that can be run individually or
via the end-to-end `gf_run_pipeline()` wrapper.

---

## Stage 0 — Initialisation

### `gf_init(envname, python, reinstall)`
Creates (or reuses) a conda environment named `"genoformer"`, installs PyTorch +
scikit-learn + numpy, and sources `model.py` / `train.py` into the Python session.
Must be called once per R session before any modelling function.

```r
gf_init()                        # uses "genoformer" conda env
gf_init(reinstall = TRUE)        # force-reinstall Python deps
```

---

## Stage 1 — Data acquisition & QC

### `gf_download_1kg(outdir, chroms, panel_only, ncores)`
Downloads per-chromosome VCFs and the sample panel file from the EBI 1KGP3 FTP.
Parallelised across chromosomes. Returns paths of downloaded files.

### `gf_qc(vcf_dir, pgs_file, outdir, maf, hwe, geno_miss, mind_miss, ld_window, ld_step, ld_r2, plink2, chroms)`
Full PLINK2 genotype QC pipeline:
1. VCF → pgen per chromosome
2. Merge chromosomes
3. Apply MAF / HWE / missingness filters
4. LD prune (for ancestry PCA)
5. Intersect with PGS scoring file variants

Returns a list with `pgen_prefix`, `pruned_prefix`, `dosage_raw`, `n_samples`,
`n_snps_qc`, `n_snps_pgs`.

---

## Stage 2 — Global ancestry inference

### `gf_ancestry(dosage_mat, panel_file, n_pcs, method, n_trees, seed)`
Infers genome-wide superpopulation ancestry for all samples.

- `method = "r"` — `prcomp` + `ranger` Random Forest (no Python needed)
- `method = "python"` — sklearn PCA + gradient boosting

Returns a `gf_ancestry` S3 object with:
- `pcs` — (N × n_pcs) PC matrix
- `anc_probs` — (N × 5) soft ancestry proportions (AFR, AMR, EAS, EUR, SAS)
- `labels` — MAP superpopulation per sample
- `entropy` — admixture entropy per sample (high = admixed)

### `print(anc)` / `plot(anc)`
Print a summary table or PCA scatter coloured by inferred superpopulation.

---

## Stage 3 — Local ancestry inference

### `gf_local_ancestry(method, ...)` — unified dispatcher
Single entry point. Dispatches to the appropriate backend and returns a
**(N × M × K) array** — one K-vector of ancestry posteriors per sample per SNP.
All downstream functions (`gf_train`, `gf_predict`) accept this array directly.

| `method` | Tool | Phasing needed | Reference panel |
|---|---|---|---|
| `"flare"` | FLARE Java jar | No | Yes |
| `"rfmix2"` | rfmix2 binary | Yes (SHAPEIT4) | Yes |
| `"windowed_pca"` | pure R | No | No |

```r
# FLARE (recommended)
la <- gf_local_ancestry("flare", vcf_file=..., ref_vcf=..., sample_map=...,
                          gmap_dir=..., snp_ids=..., snp_positions=...,
                          snp_chroms=..., sample_ids=..., global_anc=anc$anc_probs)

# Windowed PCA (no external tools)
la <- gf_local_ancestry("windowed_pca", dosage_mat=..., snp_chroms=...)
```

### `gf_local_ancestry_flare(vcf_file, ref_vcf, sample_map, gmap_dir, outdir, chroms, k_pops, flare_jar, java, nthreads, em_its)`
Runs FLARE per chromosome. Accepts unphased VCF directly. Returns named list of
output prefixes (one per chromosome). Pass to `gf_parse_flare()` to get the array.

### `gf_parse_flare(flare_prefixes, snp_ids, snp_positions, snp_chroms, sample_ids, global_anc, entropy_threshold)`
Parses FLARE `.anc.vcf.gz` files. Averages hap0 + hap1 posteriors per SNP.
Fills non-admixed samples (entropy < `entropy_threshold`) with global ancestry
as fallback. Returns **(N × M × K) array**.

### `gf_windowed_pca(dosage_mat, snp_chroms, k_pops, window_size, step_size, min_var, scale_snps, n_cores)`
Pure R sliding-window PCA proxy for local ancestry. Within each window of
`window_size` SNPs, computes top `k_pops` PCs and softmax-normalises scores to
produce pseudo-ancestry proportions summing to 1. Supports overlapping windows
(`step_size < window_size`) for smoother estimates. Returns **(N × M × k_pops) array**.

### `gf_phase(pgen_prefix, outdir, gmap_dir, chroms, shapeit4, scaffold, ncores)`
SHAPEIT4 wrapper — required pre-processing for `method = "rfmix2"` only.
Not needed for FLARE or windowed PCA.

### `gf_local_ancestry_rfmix2(vcf_file, ref_vcf, sample_map, gmap_dir, outdir, chroms, k_pops, rfmix2, em_iterations)`
Legacy RFMix2 wrapper. Requires phased VCFs. Pass output to `gf_align_rfmix2()`.

### `gf_align_rfmix2(rfmix_prefixes, snp_ids, snp_positions, snp_chroms, sample_ids, global_anc)`
Parses RFMix2 `.fb.tsv` forward-backward posteriors and aligns to dosage SNP
order via nearest-SNP interpolation. Returns **(N × M × K) array**.

### `gf_local_anc_summary(local_anc_arr, snp_chroms, pop_labels)`
Diagnostic summary: mean ancestry per labelled population per chromosome. Useful
for checking that local ancestry inference is sensible before training.

---

## Stage 4 — Model construction

### `gf_build_model(n_snps, n_ld_blocks, d_model, n_heads, n_layers, anc_dim, k_pops, d_local, use_local_anc, dropout, device)`
Constructs a `GenoFormer` PyTorch module via reticulate. Returns a `gf_model`
S3 object.

Key architecture parameters:
- `d_model` — transformer hidden dim (default 256)
- `n_layers` — transformer depth (default 6)
- `use_local_anc` — if TRUE, fuses local ancestry into each SNP token embedding
- `k_pops` — K from local ancestry array (must match `dim(local_anc)[3]`)
- `d_local` — projection dim for local ancestry embedding (default 64)

```r
# Dual-level (local + global ancestry)
model <- gf_build_model(n_snps = 80000, k_pops = 5, d_local = 64,
                         use_local_anc = TRUE)

# Global-only fallback
model <- gf_build_model(n_snps = 80000, use_local_anc = FALSE)
```

### `print(model)`
Prints device, parameter count, and full architecture config.

---

## Stage 5 — Training

### `gf_train(model, dosage_mat, anc_probs, pgs_weights, ld_blocks, chroms, phenotype, pop_labels, local_anc, epochs, batch_size, lr, weight_decay, anc_weight, calib_weight, val_split, save_path, verbose)`
Multi-task training via AdamW + cosine LR schedule. Three loss terms:

1. **PRS MSE** — primary regression loss
2. **Population classification CE** × `anc_weight` — auxiliary ancestry head
3. **Ancestry calibration loss** × `calib_weight` — penalises heteroscedastic
   residuals across superpopulations

`local_anc` is optional — if NULL the model trains in global-only mode.
Returns `list(model, history)` where `history` is a data frame of losses per epoch.

---

## Stage 6 — Inference

### `gf_predict(model, dosage_mat, anc_probs, pgs_weights, ld_blocks, chroms, effect_alleles, dosage_alleles, local_anc, batch_size)`
Batched inference. Returns a `gf_prs` tibble with:
- `transformer_prs` — standardised GenoFormer PRS
- `classical_prs` — standardised classical weighted-dosage sum (C+T)
- `local_anc_available` — logical flag
- One column per ancestry population (global soft proportions)

If the model was built with `use_local_anc=TRUE` but `local_anc=NULL` is passed,
a warning is issued and global ancestry is broadcast per-SNP as a fallback.

---

## Stage 7 — Evaluation

### `gf_evaluate(prs_tbl, phenotype, pop_labels, populations, ref_pop, n_bootstrap)`
Cross-ancestry portability evaluation following Martin et al. (2019). Computes per
population:
- Pearson R² with bootstrap 95% CIs
- Calibration slope (observed ~ predicted)
- Portability ratio: R²_pop / R²_ref_pop

Returns a `gf_eval` list with `metrics` (data frame), `portability_plot`
(ggplot2), and `distribution_plot` (PRS density by population).

---

## Model persistence

### `gf_save_model(model, path, save_config)`
Saves PyTorch state dict to `.pt` file. If `save_config=TRUE`, also writes a
`_config.json` alongside.

### `gf_load_model(path, config, device)`
Loads weights from `.pt` file. If `config=NULL`, looks for `_config.json`
next to the weights file to reconstruct architecture automatically.

---

## End-to-end pipeline

### `gf_run_pipeline(dosage_raw, pgs_file, panel_file, phenotype, outdir, local_anc_method, vcf_file, ref_vcf, sample_map, gmap_dir, k_pops, flare_jar, run_phasing, pgen_prefix, model_config, train_config, ancestry_method, device, chroms, seed)`
Runs all 6 stages in sequence. `local_anc_method` controls local ancestry:
- `"flare"` — FLARE (recommended; requires `vcf_file`, `ref_vcf`, `sample_map`, `gmap_dir`)
- `"windowed_pca"` — pure R, no external tools
- `"rfmix2"` — legacy; add `run_phasing=TRUE` if VCFs are unphased
- `"none"` — global-only conditioning

Returns a named list: `ancestry`, `local_anc_arr`, `model`, `history`, `prs`,
`evaluation`, `paths`.

---

## S3 classes produced

| Class | Created by | Key fields |
|---|---|---|
| `gf_ancestry` | `gf_ancestry()` | `pcs`, `anc_probs`, `labels`, `entropy` |
| `gf_model` | `gf_build_model()` | `model_py`, `device`, `n_params`, `config` |
| `gf_prs` | `gf_predict()` | tibble with `transformer_prs`, `classical_prs`, ancestry cols |
| `gf_eval` | `gf_evaluate()` | `metrics`, `portability_plot`, `distribution_plot` |

---

## External dependencies

| Tool | Used by | Required for |
|---|---|---|
| PLINK2 | `gf_qc()` | Genotype QC |
| conda / Miniconda | `gf_init()` | Python backend |
| PyTorch | `gf_build_model()`, `gf_train()`, `gf_predict()` | Model training & inference |
| FLARE (`flare.jar`) | `gf_local_ancestry_flare()` | Local ancestry (recommended) |
| SHAPEIT4 | `gf_phase()` | Pre-phasing for RFMix2 only |
| RFMix2 | `gf_local_ancestry_rfmix2()` | Local ancestry (legacy) |
| ranger | `gf_ancestry(method="r")` | Global ancestry (R-only path) |
