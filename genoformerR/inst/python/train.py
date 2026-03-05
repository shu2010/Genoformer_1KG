# genoformerR/inst/python/train.py  (v0.2.0)
# Training loop for GenoFormer — handles both global-only and dual-level ancestry

import torch
import torch.nn.functional as F


def prs_loss(pred, true):
    return F.mse_loss(pred, true)


def pop_loss(logits, labels):
    return F.cross_entropy(logits, labels)


def ancestry_calibration_loss(pred_prs, true_prs, anc_probs, n_pops=5):
    """
    Penalise heteroscedastic residuals across populations.
    Works with either global (B, K) or local (mean-collapsed to B, K) ancestry.
    """
    residuals = pred_prs - true_prs
    total = torch.tensor(0.0, device=pred_prs.device)
    for k in range(n_pops):
        w     = anc_probs[:, k]
        w_sum = w.sum() + 1e-8
        mu_r  = (w * residuals).sum() / w_sum
        var_r = (w * (residuals - mu_r) ** 2).sum() / w_sum
        total = total + var_r
    return total / n_pops


def _iter_batches(data_dict, batch_size, shuffle=True):
    n   = data_dict["dosage"].shape[0]
    idx = torch.randperm(n) if shuffle else torch.arange(n)
    for start in range(0, n, batch_size):
        bi = idx[start: start + batch_size]
        yield {k: v[bi] for k, v in data_dict.items()}


def train_genoformer(model, tr, va,
                     epochs=100, batch_size=64,
                     lr=1e-4, weight_decay=0.01,
                     anc_weight=0.3, calib_weight=0.5,
                     save_path=None, verbose=True):
    """
    Train GenoFormer (v0.2).

    tr / va are dicts with keys:
        dosage, log_beta, ld_block, chrom,
        global_anc,           # (B, anc_dim)
        local_anc,            # (B, L, K) or None
        prs, pop

    local_anc key is optional — if absent or None, model runs in global-only mode.
    """
    opt   = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=weight_decay)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs, eta_min=1e-6)

    has_local = "local_anc" in tr and tr["local_anc"] is not None

    best_val = float("inf")
    history  = {"epoch": [], "train_loss": [], "val_loss": [], "lr": []}

    for epoch in range(1, epochs + 1):
        # ── Train ─────────────────────────────────────────────────────────────
        model.train()
        tr_loss, n_b = 0.0, 0
        for batch in _iter_batches(tr, batch_size, shuffle=True):
            local_anc_b = batch.get("local_anc")   # None if not present

            pred_prs, pop_logits = model(
                batch["dosage"], batch["log_beta"],
                batch["ld_block"], batch["chrom"],
                batch["global_anc"], local_anc_b
            )
            # Calibration uses global ancestry for interpretability
            L = (prs_loss(pred_prs, batch["prs"])
                 + anc_weight   * pop_loss(pop_logits, batch["pop"])
                 + calib_weight * ancestry_calibration_loss(
                       pred_prs, batch["prs"], batch["global_anc"]))
            opt.zero_grad()
            L.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            tr_loss += L.item(); n_b += 1
        tr_loss /= max(n_b, 1)

        # ── Validate ──────────────────────────────────────────────────────────
        model.eval()
        va_loss, n_v = 0.0, 0
        with torch.no_grad():
            for batch in _iter_batches(va, batch_size, shuffle=False):
                local_anc_b = batch.get("local_anc")
                pred_prs, pop_logits = model(
                    batch["dosage"], batch["log_beta"],
                    batch["ld_block"], batch["chrom"],
                    batch["global_anc"], local_anc_b
                )
                L = (prs_loss(pred_prs, batch["prs"])
                     + anc_weight   * pop_loss(pop_logits, batch["pop"])
                     + calib_weight * ancestry_calibration_loss(
                           pred_prs, batch["prs"], batch["global_anc"]))
                va_loss += L.item(); n_v += 1
        va_loss /= max(n_v, 1)
        sched.step()
        cur_lr = sched.get_last_lr()[0]

        history["epoch"].append(epoch)
        history["train_loss"].append(round(tr_loss, 6))
        history["val_loss"].append(round(va_loss, 6))
        history["lr"].append(round(cur_lr, 8))

        if va_loss < best_val:
            best_val = va_loss
            if save_path:
                torch.save(model.state_dict(), save_path)

        if verbose and epoch % 10 == 0:
            mode = "dual" if has_local else "global"
            print(f"[{mode}] Epoch {epoch:4d} | "
                  f"train={tr_loss:.4f} | val={va_loss:.4f} | lr={cur_lr:.2e}")

    return {"history": history, "best_val_loss": best_val}
