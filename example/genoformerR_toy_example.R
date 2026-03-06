# =============================================================================
# genoformerR — Toy Example (pure R, no Python / PLINK2 required)
#
# Demonstrates every major R-side function using fully synthetic data.
# Python-dependent functions (gf_build_model, gf_train, gf_predict) are
# shown with working calls but wrapped in tryCatch so the script runs
# completely even without the PyTorch backend installed.
#
# Run time: ~30 seconds (R-only path), ~3 min (with Python)
# R packages needed: genoformerR, ranger, ggplot2, data.table, dplyr, tibble,
#                    cli, glue, fs
# =============================================================================

# ── 0. Install / load ─────────────────────────────────────────────────────────

# Install from the zip if not yet installed:
# install.packages("path/to/genoformerR_pkg.zip", repos=NULL, type="source")

library(genoformerR)
library(ggplot2)

set.seed(42)

# =============================================================================
# STEP 1  Simulate a toy dataset
# =============================================================================
# We simulate:
#   N = 200 samples from 5 superpopulations (AFR, AMR, EAS, EUR, SAS)
#   M = 500 SNPs spread across chromosomes 1–5

N_per_pop <- 40L          # samples per population
N         <- N_per_pop * 5L
M         <- 500L
CHROMS    <- rep(1:5, each = M / 5)        # 100 SNPs per chromosome
POP_NAMES <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# -- Population-specific allele frequencies ------------------------------------
# Each pop has slightly different minor allele frequencies
pop_mafs <- list(
  AFR = runif(M, 0.05, 0.50),
  AMR = runif(M, 0.05, 0.45),
  EAS = runif(M, 0.05, 0.40),
  EUR = runif(M, 0.05, 0.48),
  SAS = runif(M, 0.05, 0.42)
)

# -- Simulate dosage matrix (HWE random draws) --------------------------------
sim_dosage <- function(maf, n) {
  p <- maf
  sample(c(0L, 1L, 2L), n, replace = TRUE,
         prob = c((1-p)^2, 2*p*(1-p), p^2))
}

pop_labels <- rep(POP_NAMES, each = N_per_pop)
dosage_mat <- do.call(rbind, lapply(seq_along(POP_NAMES), function(pi) {
  pop   <- POP_NAMES[pi]
  mafs  <- pop_mafs[[pop]]
  do.call(rbind, lapply(seq_len(N_per_pop), function(i)
    vapply(mafs, function(maf) sim_dosage(maf, 1L), integer(1L))
  ))
}))
sample_ids <- paste0("sample_", seq_len(N))
snp_ids    <- paste0("chr", CHROMS, "_", seq(1000L, by = 1000L, length.out = M),
                     "_A_T")
rownames(dosage_mat) <- sample_ids
colnames(dosage_mat) <- snp_ids
storage.mode(dosage_mat) <- "numeric"

cat(sprintf("Dosage matrix: %d samples × %d SNPs\n", N, M))

# -- Simulate PGS effect weights (sparse causal model) ------------------------
# 50 causal SNPs with population-stratified effects
n_causal     <- 50L
causal_idx   <- sort(sample.int(M, n_causal))
beta_true    <- numeric(M)
beta_true[causal_idx] <- rnorm(n_causal, sd = 0.3)

# -- Simulate phenotype with ancestry-stratified noise ------------------------
# True PRS + pop-specific environmental effect + noise
pop_env_effects <- c(AFR = -0.5, AMR = -0.2, EAS = 0.1, EUR = 0.3, SAS = 0.0)
env_effect <- pop_env_effects[pop_labels]
prs_true   <- as.numeric(dosage_mat %*% beta_true)
phenotype  <- scale(prs_true)[, 1] + env_effect + rnorm(N, sd = 0.8)

cat(sprintf("Phenotype: mean=%.2f, sd=%.2f\n",
            mean(phenotype), sd(phenotype)))

# -- Simulate 1KGP3 panel file ------------------------------------------------
panel <- data.frame(
  sample    = sample_ids,
  pop       = pop_labels,
  super_pop = pop_labels,
  gender    = sample(c("male", "female"), N, replace = TRUE),
  stringsAsFactors = FALSE
)
panel_path <- tempfile(fileext = ".panel")
data.table::fwrite(panel, panel_path, sep = "\t")

# -- SNP metadata helpers ------------------------------------------------------
snp_positions <- seq(1000L, by = 1000L, length.out = M)
snp_chroms    <- CHROMS
ld_blocks     <- as.integer(factor(
  paste(snp_chroms, ceiling(seq_along(snp_ids) / 50L), sep = "_")
)) - 1L
pgs_weights   <- beta_true   # in practice, loaded from PGS Catalog


# =============================================================================
# STEP 2  Global ancestry inference  (pure R — no Python needed)
# =============================================================================
cat("\n── Step 2: Global ancestry inference ──\n")

anc <- gf_ancestry(
  dosage_mat = dosage_mat,
  panel_file = panel_path,
  n_pcs      = 10L,
  method     = "r",       # prcomp + ranger; no Python
  n_trees    = 100L       # reduced for speed in toy run
)

print(anc)

# Visualise: PCA scatter coloured by inferred superpopulation
pca_plot <- plot(anc)
print(pca_plot)

# Admixture entropy (higher = more admixed)
cat("Admixture entropy summary:\n")
print(summary(anc$entropy))

# Soft ancestry proportions — first few rows
cat("\nSoft ancestry probabilities (first 6 samples):\n")
print(round(head(anc$anc_probs), 3))


# =============================================================================
# STEP 3  Local ancestry — windowed PCA  (pure R — no external tools)
# =============================================================================
cat("\n── Step 3: Local ancestry (windowed PCA) ──\n")

local_anc_wpca <- gf_windowed_pca(
  dosage_mat  = dosage_mat,
  snp_chroms  = snp_chroms,
  k_pops      = 5L,           # pseudo-ancestry dimensions
  window_size = 50L,          # SNPs per window (small for toy data)
  step_size   = 25L           # 50% overlap
)

cat(sprintf("Local ancestry array dims: [%s]\n",
            paste(dim(local_anc_wpca), collapse = " × ")))

# Row-sum check: each (sample, SNP) slice should sum to 1
row_sums <- apply(local_anc_wpca, c(1, 2), sum)
cat(sprintf("Row sums: min=%.6f, max=%.6f (should be ~1)\n",
            min(row_sums), max(row_sums)))

# Summary: mean pseudo-ancestry per labelled population per chromosome
la_summary <- gf_local_anc_summary(local_anc_wpca, snp_chroms, pop_labels)
cat("\nLocal ancestry summary (first 10 rows):\n")
print(head(la_summary, 10))

# Also demonstrate the unified dispatcher
cat("\nSame result via unified dispatcher:\n")
local_anc_arr <- gf_local_ancestry(
  method      = "windowed_pca",
  dosage_mat  = dosage_mat,
  snp_chroms  = snp_chroms,
  k_pops      = 5L,
  window_size = 50L,
  step_size   = 25L
)
cat(sprintf("Dispatcher output dims: [%s]\n",
            paste(dim(local_anc_arr), collapse = " × ")))


# =============================================================================
# STEP 4  Model build, training & inference  (requires Python / PyTorch)
# =============================================================================
cat("\n── Step 4: Python backend (GenoFormer model) ──\n")

python_ok <- tryCatch({
  gf_init()
  TRUE
}, error = function(e) {
  cat("Python backend not available:", conditionMessage(e), "\n")
  cat("Steps 4–6 will be skipped. To enable, run gf_init() after installing",
      "PyTorch in a conda env named 'genoformer'.\n")
  FALSE
})

if (python_ok) {

  # ── 4a. Build model (small config for quick toy run) ──────────────────────
  model <- gf_build_model(
    n_snps        = M,
    n_ld_blocks   = max(ld_blocks) + 1L,
    d_model       = 64L,     # tiny for demo
    n_heads       = 4L,
    n_layers      = 2L,
    anc_dim       = 5L,
    k_pops        = 5L,      # matches dim(local_anc_arr)[3]
    d_local       = 16L,     # local ancestry projection dim
    use_local_anc = TRUE,
    dropout       = 0.1
  )
  print(model)

  # ── 4b. Also show global-only model for comparison ─────────────────────────
  model_global <- gf_build_model(
    n_snps      = M,
    n_ld_blocks = max(ld_blocks) + 1L,
    d_model     = 64L, n_heads = 4L, n_layers = 2L,
    use_local_anc = FALSE
  )
  cat(sprintf("Parameter count — dual-level: %s  |  global-only: %s\n",
              format(model$n_params, big.mark = ","),
              format(model_global$n_params, big.mark = ",")))

  # ── 4c. Train (very short for toy purposes) ────────────────────────────────
  pop_int <- as.integer(factor(pop_labels,
                               levels = POP_NAMES)) - 1L

  trained <- gf_train(
    model        = model,
    dosage_mat   = dosage_mat,
    anc_probs    = anc$anc_probs,
    pgs_weights  = pgs_weights,
    ld_blocks    = ld_blocks,
    chroms       = snp_chroms,
    phenotype    = phenotype,
    pop_labels   = pop_int,
    local_anc    = local_anc_arr,   # (N × M × K) from windowed PCA
    epochs       = 10L,             # low for demo; use 100+ in practice
    batch_size   = 32L,
    lr           = 1e-3,
    val_split    = 0.15,
    save_path    = tempfile(fileext = ".pt"),
    verbose      = TRUE
  )
  model   <- trained$model
  history <- trained$history

  cat("\nTraining history (last 5 epochs):\n")
  print(tail(history, 5))

  # Plot training curve
  history_plot <- ggplot(history, aes(x = epoch)) +
    geom_line(aes(y = train_loss, colour = "Train")) +
    geom_line(aes(y = val_loss,   colour = "Val")) +
    scale_colour_manual(values = c(Train = "#2E75B6", Val = "#E36C09"),
                        name = NULL) +
    labs(title = "GenoFormer training loss",
         x = "Epoch", y = "Loss") +
    theme_minimal(base_size = 12)
  print(history_plot)


  # ── 4d. Save & reload model ────────────────────────────────────────────────
  model_path <- tempfile(fileext = ".pt")
  gf_save_model(model, path = model_path, save_config = TRUE)
  cat(sprintf("Model saved to: %s\n", model_path))

  model_loaded <- gf_load_model(model_path)
  cat("Model reloaded successfully.\n")


  # ── 5. PRS inference ────────────────────────────────────────────────────────
  cat("\n── Step 5: PRS computation ──\n")

  prs <- gf_predict(
    model          = model,
    dosage_mat     = dosage_mat,
    anc_probs      = anc$anc_probs,
    pgs_weights    = pgs_weights,
    ld_blocks      = ld_blocks,
    chroms         = snp_chroms,
    local_anc      = local_anc_arr,
    batch_size     = 64L
  )
  prs$population <- pop_labels
  prs$true_pheno <- phenotype

  cat("\nPRS output (first 6 rows):\n")
  print(head(prs[, c("transformer_prs", "classical_prs",
                      "local_anc_available", "population")]))

  # Quick correlation check
  r_transformer <- cor(prs$transformer_prs, phenotype)
  r_classical   <- cor(prs$classical_prs,   phenotype)
  cat(sprintf("\nCorrelation with phenotype:\n"))
  cat(sprintf("  Transformer PRS: r = %.3f\n", r_transformer))
  cat(sprintf("  Classical PRS:   r = %.3f\n", r_classical))


  # ── 6. Portability evaluation ───────────────────────────────────────────────
  cat("\n── Step 6: Portability evaluation ──\n")

  eval_res <- gf_evaluate(
    prs_tbl    = prs,
    phenotype  = phenotype,
    pop_labels = pop_labels,
    ref_pop    = "EUR",
    n_bootstrap = 200L   # low for speed
  )

  cat("\nPortability metrics:\n")
  print(eval_res$metrics)

  # Show the two evaluation plots
  print(eval_res$portability_plot)
  print(eval_res$distribution_plot)

}  # end python_ok block


# =============================================================================
# STEP 7  Classical PRS baseline  (pure R, no Python)
# =============================================================================
cat("\n── Step 7: Classical PRS baseline (pure R) ──\n")

# Weighted dosage sum — C+T equivalent
classical_prs <- as.numeric(dosage_mat %*% pgs_weights)
classical_prs <- scale(classical_prs)[, 1]

# R² per population
cat("Classical PRS R² per population:\n")
for (pop in POP_NAMES) {
  mask <- pop_labels == pop
  r2   <- cor(classical_prs[mask], phenotype[mask])^2
  cat(sprintf("  %-4s: R² = %.4f\n", pop, r2))
}


# =============================================================================
# STEP 8  Demonstrate individual stage functions in isolation
# =============================================================================
cat("\n── Step 8: Individual function demos ──\n")

# 8a. Windowed PCA with overlapping windows
cat("gf_windowed_pca with step_size < window_size (overlapping):\n")
la_overlap <- gf_windowed_pca(dosage_mat, snp_chroms,
                               k_pops = 3, window_size = 40, step_size = 20)
cat(sprintf("  dims: [%s]\n", paste(dim(la_overlap), collapse=" × ")))

# 8b. Local ancestry summary across chromosomes
cat("\ngf_local_anc_summary:\n")
summ <- gf_local_anc_summary(local_anc_arr, snp_chroms, pop_labels)
print(head(summ, 6))

# 8c. Show that the dispatcher rejects missing args for FLARE
cat("\ngf_local_ancestry dispatcher error handling:\n")
tryCatch(
  gf_local_ancestry(method = "flare", dosage_mat = dosage_mat),
  error = function(e) cat("Expected error caught:", conditionMessage(e), "\n")
)

# 8d. Confirm internal helpers work as expected
cat("\nInternal SNP metadata parsers:\n")
test_ids <- c("chr1_1000_A_T", "chr3_5000_G_C", "chr5_99000_T_A")
cat("  chromosomes: ", paste(genoformerR:::.snp_chroms(test_ids), collapse=" "), "\n")
cat("  positions:   ", paste(genoformerR:::.snp_positions(test_ids), collapse=" "), "\n")
blocks_test <- genoformerR:::.assign_ld_blocks(test_ids,
                                               genoformerR:::.snp_chroms(test_ids))
cat("  LD blocks:   ", paste(blocks_test, collapse=" "), "\n")


# =============================================================================
# SESSION INFO
# =============================================================================
cat("\n── Session info ──\n")
cat(sprintf("R version: %s\n", R.version$version.string))
cat(sprintf("genoformerR version: %s\n",
            as.character(packageVersion("genoformerR"))))
if (python_ok) {
  torch_ver <- tryCatch(
    reticulate::py_to_r(reticulate::import("torch")$`__version__`),
    error = function(e) "unavailable"
  )
  cat(sprintf("PyTorch version: %s\n", torch_ver))
}
cat("\nToy example complete.\n")
