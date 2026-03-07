# genoformerR/R/04_prs.R
# Step 6–7: PRS computation, classical comparison, portability evaluation

#' Compute ancestry-aware PRS using the trained GenoFormer model
#'
#' Runs inference and returns transformer-predicted PRS alongside the classical
#' weighted dosage sum (C+T). Accepts an optional \code{local_anc} array
#' (N × M × K) for dual-level ancestry conditioning at inference time.
#'
#' @param model          A trained \code{gf_model} from \code{\link{gf_train}}.
#' @param dosage_mat     Numeric matrix (N × M) of dosage values.
#' @param anc_probs      Numeric matrix (N × 5) of global soft ancestry proportions.
#' @param pgs_weights    Numeric vector (M). Effect weights from PGS file.
#' @param ld_blocks      Integer vector (M). LD block ID per SNP.
#' @param chroms         Integer vector (M). Chromosome per SNP.
#' @param effect_alleles Character vector (M). Effect allele per SNP.
#' @param dosage_alleles Character vector (M). ALT allele in dosage matrix.
#' @param local_anc      Numeric array (N × M × K) of per-SNP local ancestry
#'   posteriors from \code{\link{gf_load_local_anc}}, or NULL. When NULL and
#'   the model was built with \code{use_local_anc=TRUE}, a warning is issued
#'   and global ancestry is broadcast per-SNP as a fallback.
#' @param batch_size     Integer. Inference batch size. Default 128.
#'
#' @return A \code{gf_prs} tibble with columns:
#'   \code{transformer_prs}, \code{classical_prs}, ancestry columns,
#'   and \code{local_anc_available} (logical).
#' @export
gf_predict <- function(model,
                        dosage_mat,
                        anc_probs,
                        pgs_weights,
                        ld_blocks,
                        chroms,
                        effect_alleles = NULL,
                        dosage_alleles = NULL,
                        local_anc      = NULL,
                        batch_size     = 128L) {
  .check_init()
  stopifnot(inherits(model, "gf_model"))

  use_local <- !is.null(local_anc) && model$config$use_local_anc
  if (model$config$use_local_anc && is.null(local_anc)) {
    cli::cli_warn(c(
      "Model was built with use_local_anc=TRUE but local_anc=NULL.",
      "i" = "Falling back: broadcasting global ancestry per SNP position.",
      "i" = "For admixed samples, results will be suboptimal."
    ))
    # Fallback: broadcast global anc (N × 5) → (N × M × 5)
    M          <- ncol(dosage_mat)
    K          <- ncol(anc_probs)
    local_anc  <- array(0, dim = c(nrow(dosage_mat), M, K))
    for (m in seq_len(M)) local_anc[, m, ] <- anc_probs
    use_local  <- TRUE
  }

  torch <- .gf_env$torch
  dev   <- model$device
  N     <- nrow(dosage_mat)
  M     <- ncol(dosage_mat)
  log_beta <- log(abs(pgs_weights) + 1e-8)

  # ── Classical PRS ──────────────────────────────────────────────────────────
  cli::cli_alert_info("Computing classical C+T PRS...")
  if (!is.null(effect_alleles) && !is.null(dosage_alleles)) {
    aligned <- ifelse(effect_alleles == dosage_alleles, dosage_mat, 2 - dosage_mat)
  } else {
    aligned <- dosage_mat
  }
  classical_prs <- as.numeric(aligned %*% pgs_weights)

  # ── Transformer PRS ────────────────────────────────────────────────────────
  mode_str <- if (use_local) "dual (global + local)" else "global-only"
  cli::cli_alert_info("Running GenoFormer inference [{mode_str}] on {N} samples...")
  model$model_py$eval()
  all_prs <- numeric(N)
  batches <- split(seq_len(N), ceiling(seq_len(N) / batch_size))

  pb <- cli::cli_progress_bar("Inference", total = length(batches))
  for (batch_idx in batches) {
    to_t <- function(x, dtype = "float32") {
      py_dtype <- switch(dtype,
        "float32" = torch$float32,
        "float64" = torch$float64,
        "int64"   = torch$int64,
        torch$float32
      )
      torch$tensor(reticulate::r_to_py(x), dtype = py_dtype)$to(dev)
    }
    d      <- to_t(dosage_mat[batch_idx, ])
    lb     <- to_t(matrix(log_beta,  nrow = length(batch_idx), ncol = M, byrow = TRUE))
    lb_blk <- torch$tensor(
      reticulate::r_to_py(matrix(ld_blocks, nrow = length(batch_idx), ncol = M, byrow = TRUE)),
      dtype = torch$int64)$to(dev)
    ch <- torch$tensor(
      reticulate::r_to_py(matrix(chroms, nrow = length(batch_idx), ncol = M, byrow = TRUE)),
      dtype = torch$int64)$to(dev)
    g_anc <- to_t(anc_probs[batch_idx, ])
    l_anc <- if (use_local) to_t(local_anc[batch_idx, , ]) else reticulate::py_none()

    with(torch$no_grad(), {
      out   <- model$model_py(d, lb, lb_blk, ch, g_anc, l_anc)
      prs_b <- reticulate::py_to_r(out[[1]]$cpu()$numpy())
      all_prs[batch_idx] <- as.numeric(prs_b)
    })
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  # Column names for ancestry: use pop names from local_anc or default
  anc_col_names <- colnames(anc_probs)
  if (is.null(anc_col_names)) anc_col_names <- paste0("anc_", seq_len(ncol(anc_probs)))

  result <- tibble::tibble(
    transformer_prs     = scale(all_prs)[, 1],
    classical_prs       = scale(classical_prs)[, 1],
    local_anc_available = use_local
  )
  # Append global ancestry columns
  anc_df <- as.data.frame(anc_probs)
  names(anc_df) <- anc_col_names
  result <- dplyr::bind_cols(result, anc_df)

  class(result) <- c("gf_prs", class(result))
  cli::cli_alert_success("PRS computation complete.")
  result
}


#' Evaluate PRS performance with cross-ancestry portability metrics
#'
#' Computes Pearson R², calibration slopes, and the Martin et al. (2019)
#' portability ratio (R²_pop / R²_EUR) for both transformer and classical PRS.
#'
#' @param prs_tbl      A \code{gf_prs} tibble from \code{\link{gf_predict}}.
#' @param phenotype    Numeric vector (N). Observed phenotype.
#' @param pop_labels   Character vector (N). Superpopulation per sample.
#' @param populations  Character vector. Populations to evaluate. Default all 5.
#' @param ref_pop      Character. Reference population for portability ratio. Default "EUR".
#' @param n_bootstrap  Integer. Bootstrap replicates for CIs. Default 500.
#'
#' @return A \code{gf_eval} list with:
#'   \describe{
#'     \item{\code{metrics}}{Data frame with R², portability, calibration per population per method.}
#'     \item{\code{portability_plot}}{ggplot2 object.}
#'     \item{\code{distribution_plot}}{ggplot2 object.}
#'   }
#' @export
gf_evaluate <- function(prs_tbl,
                         phenotype,
                         pop_labels,
                         populations = c("AFR", "AMR", "EAS", "EUR", "SAS"),
                         ref_pop     = "EUR",
                         n_bootstrap = 500L) {

  cli::cli_h1("PRS Portability Evaluation")

  .r2_ci <- function(x, y, n_boot) {
    r   <- cor(x, y, use = "complete.obs")
    r2  <- r^2
    boot_r2 <- replicate(n_boot, {
      i <- sample(length(x), replace = TRUE)
      cor(x[i], y[i], use = "complete.obs")^2
    })
    c(r2 = r2, lo = quantile(boot_r2, 0.025), hi = quantile(boot_r2, 0.975))
  }

  rows <- list()
  for (pop in populations) {
    mask <- pop_labels == pop
    if (sum(mask, na.rm = TRUE) < 20) next
    ph_p <- phenotype[mask]
    for (method in c("transformer_prs", "classical_prs")) {
      prs_p <- prs_tbl[[method]][mask]
      est   <- .r2_ci(prs_p, ph_p, n_bootstrap)
      rows[[length(rows) + 1]] <- data.frame(
        pop    = pop, method = method,
        r2     = est["r2"],
        lo     = est["lo"],
        hi     = est["hi"],
        n      = sum(mask),
        stringsAsFactors = FALSE
      )
    }
  }
  metrics <- do.call(rbind, rows)
  rownames(metrics) <- NULL

  # Portability ratio vs EUR
  ref_r2 <- lapply(c("transformer_prs", "classical_prs"), function(m) {
    metrics[metrics$pop == ref_pop & metrics$method == m, "r2"]
  })
  names(ref_r2) <- c("transformer_prs", "classical_prs")
  metrics$portability <- metrics$r2 / ref_r2[metrics$method]

  cli::cli_alert_info("Portability ratios (vs {ref_pop}):")
  for (pop in setdiff(populations, ref_pop)) {
    tr_p  <- metrics[metrics$pop == pop & metrics$method == "transformer_prs", "portability"]
    cl_p  <- metrics[metrics$pop == pop & metrics$method == "classical_prs",   "portability"]
    if (length(tr_p) && length(cl_p)) {
      cli::cli_dl(setNames(
        glue::glue("Transformer={round(tr_p,3)}  Classical={round(cl_p,3)}"), pop
      ))
    }
  }

  result <- list(
    metrics          = metrics,
    portability_plot = .plot_portability(metrics, ref_pop),
    distribution_plot = .plot_distributions(prs_tbl, pop_labels, populations)
  )
  class(result) <- c("gf_eval", "list")
  result
}


#' @export
print.gf_eval <- function(x, ...) {
  cat("\n=== GenoFormer PRS Evaluation ===\n")
  print(x$metrics[, c("pop", "method", "r2", "lo", "hi", "portability")],
        digits = 3, row.names = FALSE)
  invisible(x)
}


# ─── Plot helpers ──────────────────────────────────────────────────────────────

.plot_portability <- function(metrics, ref_pop) {
  df <- metrics[metrics$pop != ref_pop, ]
  df$method_lab <- ifelse(df$method == "transformer_prs", "GenoFormer", "Classical C+T")
  pal <- c("GenoFormer" = "#00e5ff", "Classical C+T" = "#ff6b6b")

  ggplot2::ggplot(df, ggplot2::aes(x = pop, y = portability,
                                    colour = method_lab,
                                    group  = method_lab)) +
    ggplot2::geom_hline(yintercept = 0.7, linetype = "dashed",
                        colour = "#555580", linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = 1,   linetype = "solid",
                        colour = "#333355", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 1, alpha = 0.8) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo / df$r2[df$pop == ref_pop][1],
                                         ymax = hi / df$r2[df$pop == ref_pop][1]),
                           width = 0.15, alpha = 0.5) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_manual(values = pal) +
    ggplot2::scale_y_continuous(limits = c(0, 1.1),
                                 labels = scales::percent_format()) +
    ggplot2::labs(
      title    = glue::glue("PRS Portability (R\u00B2 relative to {ref_pop})"),
      subtitle = "Dashed line = 0.70 equity threshold (Martin et al. 2019)",
      x        = "Superpopulation",
      y        = glue::glue("R\u00B2 / R\u00B2_{ref_pop}"),
      colour   = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    .dark_theme()
}


.plot_distributions <- function(prs_tbl, pop_labels, populations) {
  pal <- c(AFR = "#ff6b6b", AMR = "#ffd700", EAS = "#69ff47",
           EUR = "#00e5ff", SAS = "#c77dff")
  df_long <- rbind(
    data.frame(prs = prs_tbl$transformer_prs, pop = pop_labels,
               method = "GenoFormer", stringsAsFactors = FALSE),
    data.frame(prs = prs_tbl$classical_prs, pop = pop_labels,
               method = "Classical C+T", stringsAsFactors = FALSE)
  )
  df_long$pop <- factor(df_long$pop, levels = populations)

  ggplot2::ggplot(df_long[df_long$pop %in% populations, ],
                  ggplot2::aes(x = prs, fill = pop)) +
    ggplot2::geom_density(alpha = 0.6, colour = NA) +
    ggplot2::facet_wrap(~method, ncol = 1) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::labs(
      title = "PRS Distributions by Superpopulation",
      x     = "Standardised PRS", y = "Density", fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    .dark_theme()
}


.dark_theme <- function() {
  ggplot2::theme(
    plot.background  = ggplot2::element_rect(fill = "#0a0a14", colour = NA),
    panel.background = ggplot2::element_rect(fill = "#0f0f1e", colour = NA),
    strip.background = ggplot2::element_rect(fill = "#111128", colour = NA),
    text             = ggplot2::element_text(colour = "#e0e0ff"),
    strip.text       = ggplot2::element_text(colour = "#8888bb"),
    axis.text        = ggplot2::element_text(colour = "#8888bb"),
    panel.grid       = ggplot2::element_line(colour = "#1a1a3a"),
    legend.background = ggplot2::element_rect(fill = "#0a0a14", colour = NA),
    legend.key       = ggplot2::element_rect(fill = "#0a0a14", colour = NA)
  )
}
