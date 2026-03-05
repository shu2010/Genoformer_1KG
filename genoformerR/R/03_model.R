# genoformerR/R/03_model.R
# Step 4–5: GenoFormer model construction and training via reticulate

#' Build the GenoFormer transformer model
#'
#' Constructs a GenoFormer v0.2 instance in Python/PyTorch via reticulate.
#' Supports dual-level ancestry conditioning: global ancestry (FiLM per block)
#' and local ancestry (per-SNP token embedding from RFMix2 posteriors).
#'
#' @param n_snps        Integer. Number of PGS SNPs (sequence length). Default 80000.
#' @param n_ld_blocks   Integer. Number of LD blocks (genome-wide ~1700). Default 1700.
#' @param d_model       Integer. Transformer hidden dimension. Default 256.
#' @param n_heads       Integer. Number of attention heads. Default 8.
#' @param n_layers      Integer. Transformer depth. Default 6.
#' @param anc_dim       Integer. Global ancestry dimension (5 superpops). Default 5.
#' @param k_pops        Integer. Number of reference populations in local ancestry
#'   posteriors (K in RFMix2 output). Must match \code{dim(local_anc_arr)[3]}.
#'   Default 5.
#' @param d_local       Integer. Projection dimension for local ancestry embedding.
#'   Larger values allow richer local ancestry representations. Default 64.
#' @param use_local_anc Logical. If FALSE, disables local ancestry embedding and
#'   the model behaves identically to v0.1 (global-only). Default TRUE.
#' @param dropout       Numeric. Dropout rate. Default 0.1.
#' @param device        Character. \code{"cuda"}, \code{"mps"}, or \code{"cpu"}.
#'   Auto-detected if \code{"auto"}.
#'
#' @return An S3 object of class \code{gf_model} wrapping the PyTorch module.
#' @export
#'
#' @examples
#' \dontrun{
#' gf_init()
#' # With local ancestry (RFMix2 K=3 pops: AFR, EUR, NAT)
#' model <- gf_build_model(n_snps = 5000, k_pops = 3, d_local = 32)
#' # Global-only fallback (v0.1 behaviour)
#' model <- gf_build_model(n_snps = 5000, use_local_anc = FALSE)
#' }
gf_build_model <- function(n_snps        = 80000L,
                            n_ld_blocks   = 1700L,
                            d_model       = 256L,
                            n_heads       = 8L,
                            n_layers      = 6L,
                            anc_dim       = 5L,
                            k_pops        = 5L,
                            d_local       = 64L,
                            use_local_anc = TRUE,
                            dropout       = 0.1,
                            device        = "auto") {
  .check_init()
  device <- .resolve_device(device)
  cli::cli_alert_info("Building GenoFormer v0.2 on {.val {device}} ...")
  if (use_local_anc)
    cli::cli_alert_info("Local ancestry: k_pops={k_pops}, d_local={d_local}")
  else
    cli::cli_alert_info("Local ancestry: disabled (global-only mode)")

  torch <- .gf_env$torch
  main  <- reticulate::import_main()

  model_py <- main$GenoFormer(
    n_snps        = as.integer(n_snps),
    n_ld_blocks   = as.integer(n_ld_blocks),
    d_model       = as.integer(d_model),
    n_heads       = as.integer(n_heads),
    n_layers      = as.integer(n_layers),
    anc_dim       = as.integer(anc_dim),
    k_pops        = as.integer(k_pops),
    d_local       = as.integer(d_local),
    dropout       = dropout,
    use_local_anc = use_local_anc
  )
  model_py <- model_py$to(device)

  n_params <- reticulate::py_to_r(
    torch$tensor(
      reticulate::py_eval(
        "sum(p.numel() for p in model_py.parameters() if p.requires_grad)",
        local = list(model_py = model_py))
    )$item()
  )

  result <- structure(
    list(
      model_py  = model_py,
      device    = device,
      n_params  = n_params,
      config    = list(n_snps = n_snps, n_ld_blocks = n_ld_blocks,
                       d_model = d_model, n_heads = n_heads,
                       n_layers = n_layers, anc_dim = anc_dim,
                       k_pops = k_pops, d_local = d_local,
                       use_local_anc = use_local_anc, dropout = dropout)
    ),
    class = c("gf_model", "list")
  )
  cli::cli_alert_success(
    "GenoFormer ready: {.val {format(n_params, big.mark=',')}} parameters on {.val {device}}"
  )
  result
}


#' @export
print.gf_model <- function(x, ...) {
  cli::cli_h1("GenoFormer Model (v0.2)")
  cfg <- x$config
  cli::cli_dl(c(
    "Device"          = x$device,
    "Parameters"      = format(x$n_params, big.mark = ","),
    "SNPs (L)"        = format(cfg$n_snps,      big.mark = ","),
    "d_model"         = cfg$d_model,
    "Layers"          = cfg$n_layers,
    "Heads"           = cfg$n_heads,
    "LD blocks"       = cfg$n_ld_blocks,
    "Global anc dim"  = cfg$anc_dim,
    "Local anc (K)"   = if (cfg$use_local_anc) cfg$k_pops  else "disabled",
    "Local anc (dim)" = if (cfg$use_local_anc) cfg$d_local else "disabled"
  ))
  invisible(x)
}


#' Train the GenoFormer model
#'
#' Runs multi-task training: PRS regression (MSE) + ancestry classification
#' (cross-entropy) + ancestry calibration loss. Accepts an optional
#' \code{local_anc} array (N × M × K) from RFMix2 for dual-level conditioning.
#' When provided, local ancestry is passed as a per-SNP tensor alongside the
#' global ancestry vector.
#'
#' @param model        A \code{gf_model} object from \code{\link{gf_build_model}}.
#' @param dosage_mat   Numeric matrix (N × M) of dosage values.
#' @param anc_probs    Numeric matrix (N × 5) of global soft ancestry proportions
#'   from \code{\link{gf_ancestry}}.
#' @param pgs_weights  Numeric vector (M). Effect weights from PGS scoring file.
#' @param ld_blocks    Integer vector (M). LD block assignment per SNP.
#' @param chroms       Integer vector (M). Chromosome per SNP.
#' @param phenotype    Numeric vector (N). Continuous or binary phenotype.
#' @param pop_labels   Integer vector (N). Superpopulation integer codes (0–4).
#' @param local_anc    Numeric array (N × M × K) of per-SNP local ancestry
#'   posteriors from \code{\link{gf_load_local_anc}}, or NULL to disable local
#'   ancestry conditioning (global-only mode). Default NULL.
#' @param epochs       Integer. Training epochs. Default 100.
#' @param batch_size   Integer. Mini-batch size. Default 64.
#' @param lr           Numeric. Initial learning rate. Default 1e-4.
#' @param weight_decay Numeric. AdamW weight decay. Default 0.01.
#' @param anc_weight   Numeric. Ancestry classification loss weight. Default 0.3.
#' @param calib_weight Numeric. Calibration loss weight. Default 0.5.
#' @param val_split    Numeric. Fraction held out for validation. Default 0.1.
#' @param save_path    Character or NULL. Where to save the best checkpoint.
#' @param verbose      Logical. Print epoch-level progress. Default TRUE.
#'
#' @return A list with \code{model} and \code{history} (data frame).
#' @export
gf_train <- function(model,
                     dosage_mat,
                     anc_probs,
                     pgs_weights,
                     ld_blocks,
                     chroms,
                     phenotype,
                     pop_labels,
                     local_anc    = NULL,
                     epochs       = 100L,
                     batch_size   = 64L,
                     lr           = 1e-4,
                     weight_decay = 0.01,
                     anc_weight   = 0.3,
                     calib_weight = 0.5,
                     val_split    = 0.1,
                     save_path    = "genoformer_best.pt",
                     verbose      = TRUE) {
  .check_init()
  stopifnot(inherits(model, "gf_model"))

  use_local <- !is.null(local_anc) && model$config$use_local_anc
  if (use_local) {
    stopifnot(
      is.array(local_anc),
      length(dim(local_anc)) == 3L,
      dim(local_anc)[1] == nrow(dosage_mat),
      dim(local_anc)[2] == ncol(dosage_mat)
    )
    K_local <- dim(local_anc)[3]
    if (K_local != model$config$k_pops) {
      cli::cli_abort(
        "local_anc has K={K_local} populations but model was built with k_pops={model$config$k_pops}."
      )
    }
    cli::cli_alert_info("Dual-level conditioning: global ancestry + local ancestry (K={K_local})")
  } else {
    cli::cli_alert_info("Global-only conditioning (local_anc not provided or disabled in model)")
  }

  torch <- .gf_env$torch
  main  <- reticulate::import_main()
  dev   <- model$device

  cli::cli_h1("GenoFormer Training")
  cli::cli_alert_info("Epochs={epochs} | batch={batch_size} | lr={lr}")

  # ── Train / val split ─────────────────────────────────────────────────────
  N     <- nrow(dosage_mat)
  val_n <- max(1L, floor(N * val_split))
  idx   <- sample.int(N)
  val_i <- idx[seq_len(val_n)]
  tr_i  <- idx[-seq_len(val_n)]

  to_tensor <- function(x, dtype = "float32") {
    torch$tensor(reticulate::r_to_py(x), dtype = torch[[dtype]])$to(dev)
  }

  log_beta <- log(abs(pgs_weights) + 1e-8)

  .make_split <- function(i) {
    d <- list(
      dosage     = to_tensor(dosage_mat[i, ]),
      log_beta   = to_tensor(matrix(log_beta,   nrow = length(i), ncol = length(log_beta),   byrow = TRUE)),
      ld_block   = to_tensor(matrix(ld_blocks,  nrow = length(i), ncol = length(ld_blocks),  byrow = TRUE), "int64"),
      chrom      = to_tensor(matrix(chroms,     nrow = length(i), ncol = length(chroms),     byrow = TRUE), "int64"),
      global_anc = to_tensor(anc_probs[i, ]),
      prs        = to_tensor(scale(phenotype)[i]),
      pop        = to_tensor(pop_labels[i], "int64")
    )
    if (use_local) {
      # local_anc[i, , ] is (length(i), M, K)
      d$local_anc <- to_tensor(local_anc[i, , ])
    }
    d
  }

  tr <- .make_split(tr_i)
  va <- .make_split(val_i)

  # ── Delegate to Python training loop ──────────────────────────────────────
  result_py <- main$train_genoformer(
    model        = model$model_py,
    tr           = reticulate::r_to_py(tr),
    va           = reticulate::r_to_py(va),
    epochs       = as.integer(epochs),
    batch_size   = as.integer(batch_size),
    lr           = lr,
    weight_decay = weight_decay,
    anc_weight   = anc_weight,
    calib_weight = calib_weight,
    save_path    = if (is.null(save_path)) reticulate::py_none() else save_path,
    verbose      = verbose
  )

  history <- as.data.frame(reticulate::py_to_r(result_py$history))

  if (!is.null(save_path) && file.exists(save_path)) {
    model$model_py$load_state_dict(torch$load(save_path))
    cli::cli_alert_success("Best checkpoint loaded from {.path {save_path}}")
  }

  cli::cli_alert_success(
    "Training complete. Best val loss: {.val {round(min(history$val_loss), 4)}}"
  )
  list(model = model, history = history)
}


#' Save a trained GenoFormer model
#'
#' @param model      A \code{gf_model} object.
#' @param path       Character. Output \code{.pt} file path.
#' @param save_config Logical. Also save JSON config alongside the weights.
#' @export
gf_save_model <- function(model, path = "genoformer.pt", save_config = TRUE) {
  .check_init()
  .gf_env$torch$save(model$model_py$state_dict(), path)
  if (save_config) {
    cfg_path <- sub("\\.pt$", "_config.json", path)
    jsonlite::write_json(model$config, cfg_path, pretty = TRUE)
  }
  cli::cli_alert_success("Model saved to {.path {path}}")
  invisible(path)
}


#' Load a saved GenoFormer model
#'
#' @param path   Character. Path to \code{.pt} weights file.
#' @param config Named list or path to JSON config. If NULL, uses defaults.
#' @param device Character. Device for inference.
#' @export
gf_load_model <- function(path, config = NULL, device = "auto") {
  .check_init()
  if (is.null(config)) {
    cfg_path <- sub("\\.pt$", "_config.json", path)
    if (file.exists(cfg_path)) {
      config <- jsonlite::read_json(cfg_path)
    }
  }
  args <- if (is.null(config)) list() else config
  model <- do.call(gf_build_model, c(args, list(device = device)))
  model$model_py$load_state_dict(.gf_env$torch$load(path))
  cli::cli_alert_success("Weights loaded from {.path {path}}")
  model
}


# ─── Internal helpers ──────────────────────────────────────────────────────────
.resolve_device <- function(device) {
  if (device != "auto") return(device)
  torch <- .gf_env$torch
  if (reticulate::py_to_r(torch$cuda$is_available())) return("cuda")
  if (reticulate::py_to_r(torch$backends$mps$is_available())) return("mps")
  "cpu"
}
