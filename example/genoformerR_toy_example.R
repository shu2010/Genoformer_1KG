# =============================================================================
# genoformerR — Toy Example (v0.3.1)
#
# Demonstrates every major R-side function using fully synthetic data.
# Python-dependent functions (gf_build_model, gf_train, gf_predict) are
# shown with working calls but wrapped in tryCatch so the script runs
# completely even without the PyTorch backend.
#
# v0.3.1 changes demonstrated:
#   - gf_ancestry()     new ref_dosage_raw= parameter (Mode A vs Mode B)
#   - gf_windowed_pca() new ref_dosage_mat= and ref_window_pcas= parameters
#   - gf_save_window_pcas() / gf_load_window_pcas() new functions
#
# Run time: ~30 seconds (R-only), ~3 min (with Python)
# Packages: genoformerR, ranger, ggplot2, data.table, dplyr, tibble, cli, fs
# =============================================================================

library(genoformerR)
library(ggplot2)

set.seed(42)


# =============================================================================
# STEP 1  Simulate a toy dataset
# =============================================================================
# N = 200 samples from 5 superpopulations (40 per pop)
# M = 500 SNPs across chromosomes 1–5
# We simulate two independent cohorts:
#   ref_dosage  — stands in for 1KGP3 (labelled reference samples)
#   study_dosage — stands in for your external study cohort (unlabelled)
# This structure mirrors the real use-case for Mode A ancestry inference.

N_per_pop  <- 40L
N          <- N_per_pop * 5L
M          <- 500L
CHROMS     <- rep(1:5, each = M / 5)
POP_NAMES  <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# Population-specific allele frequencies
pop_mafs <- list(
  AFR = runif(M, 0.05, 0.50),
  AMR = runif(M, 0.05, 0.45),
  EAS = runif(M, 0.05, 0.40),
  EUR = runif(M, 0.05, 0.48),
  SAS = runif(M, 0.05, 0.42)
)

sim_dosage_mat <- function(ids_prefix, n_per_pop) {
  pop_labels_local <- rep(POP_NAMES, each = n_per_pop)
  N_local <- n_per_pop * 5L
  mat <- do.call(rbind, lapply(seq_along(POP_NAMES), function(pi) {
    pop  <- POP_NAMES[pi]
    mafs <- pop_mafs[[pop]]
    do.call(rbind, lapply(seq_len(n_per_pop), function(i)
      vapply(mafs, function(maf) {
        p <- maf
        sample(c(0L, 1L, 2L), 1L, prob = c((1-p)^2, 2*p*(1-p), p^2))
      }, integer(1L))
    ))
  }))
  rownames(mat) <- paste0(ids_prefix, "_", seq_len(N_local))
  mat
}

snp_ids    <- paste0("chr", CHROMS, "_", seq(1000L, by = 1000L, length.out = M), "_A_T")

# Reference cohort (labelled — acts as 1KGP3 stand-in)
ref_dosage        <- sim_dosage_mat("ref", N_per_pop)
colnames(ref_dosage) <- snp_ids
storage.mode(ref_dosage) <- "numeric"
ref_pop_labels    <- rep(POP_NAMES, each = N_per_pop)

# Study cohort (unlabelled — acts as your research cohort)
study_dosage      <- sim_dosage_mat("study", N_per_pop)
colnames(study_dosage) <- snp_ids
storage.mode(study_dosage) <- "numeric"
study_pop_labels  <- rep(POP_NAMES, each = N_per_pop)  # known only for eval

cat(sprintf("Reference dosage: %d x %d\n", nrow(ref_dosage), ncol(ref_dosage)))
cat(sprintf("Study dosage:     %d x %d\n", nrow(study_dosage), ncol(study_dosage)))

# Panel file — covers reference samples only (not study samples)
panel <- data.frame(
  sample    = rownames(ref_dosage),
  pop       = ref_pop_labels,
  super_pop = ref_pop_labels,
  gender    = sample(c("male","female"), nrow(ref_dosage), replace=TRUE),
  stringsAsFactors = FALSE
)
panel_path <- tempfile(fileext = ".panel")
data.table::fwrite(panel, panel_path, sep = "\t")

# SNP metadata
snp_positions <- seq(1000L, by = 1000L, length.out = M)
snp_chroms    <- CHROMS
ld_blocks     <- as.integer(factor(
  paste(snp_chroms, ceiling(seq_along(snp_ids) / 50L), sep = "_")
)) - 1L

# PGS weights (50 causal SNPs, sparse)
n_causal    <- 50L
causal_idx  <- sort(sample.int(M, n_causal))
pgs_weights <- numeric(M)
pgs_weights[causal_idx] <- rnorm(n_causal, sd = 0.3)

# Phenotype for study cohort
pop_env_effects <- c(AFR=-0.5, AMR=-0.2, EAS=0.1, EUR=0.3, SAS=0.0)
prs_true        <- as.numeric(study_dosage %*% pgs_weights)
phenotype       <- scale(prs_true)[, 1] +
                   pop_env_effects[study_pop_labels] +
                   rnorm(nrow(study_dosage), sd = 0.8)

cat(sprintf("Phenotype: mean=%.2f  sd=%.2f\n", mean(phenotype), sd(phenotype)))


# =============================================================================
# STEP 2  Global ancestry — Mode A (reference projection)
# =============================================================================
# The study cohort has no entries in the panel file (different sample IDs).
# We must provide ref_dosage as ref_dosage_raw so that:
#   - PCA is fitted on the reference samples
#   - study samples are projected into that reference PC space
#   - the Random Forest is trained on reference labels
#   - ancestry probabilities for study samples are predicted
#
# In a real run, ref_dosage_raw= would be the path to the 1KGP3 .raw file
# produced by gf_qc().  Here we pass the matrix directly.

cat("\n── Step 2: Global ancestry — Mode A ──\n")

# Write reference dosage to a temp .raw file to demonstrate the file-path API
# (In practice this file already exists from gf_qc())
ref_raw_path <- tempfile(fileext = ".raw")
ref_raw      <- cbind(
  data.frame(FID=rownames(ref_dosage), IID=rownames(ref_dosage),
             PAT=0, MAT=0, SEX=0, PHENOTYPE=-9),
  as.data.frame(ref_dosage)
)
data.table::fwrite(ref_raw, ref_raw_path, sep=" ")

anc <- gf_ancestry(
  dosage_mat     = study_dosage,     # study cohort — not in panel file
  panel_file     = panel_path,
  ref_dosage_raw = ref_raw_path,     # 1KGP3 reference dosage
  n_pcs          = 10L,
  method         = "r",
  n_trees        = 100L
)

print(anc)
pca_plot <- plot(anc)
print(pca_plot)

cat("Soft ancestry probabilities (first 6 study samples):\n")
print(round(head(anc$anc_probs), 3))

# Save — pca_obj and clf_obj inside this object are needed for scoring
anc_rds_path <- tempfile(fileext = ".rds")
saveRDS(anc, anc_rds_path)
cat(sprintf("ancestry.rds saved to: %s\n", anc_rds_path))

# Also demonstrate Mode B — valid when dosage_mat IS the reference cohort
cat("\n── Mode B (study-only, valid when dosage_mat IS reference cohort) ──\n")
anc_modeB <- gf_ancestry(
  dosage_mat = ref_dosage,           # reference samples — IDs match panel
  panel_file = panel_path
)
cat(sprintf("Mode B: mean entropy = %.3f\n", mean(anc_modeB$entropy)))

# Integer-coded population labels for training
pop_int <- as.integer(factor(study_pop_labels, levels = POP_NAMES)) - 1L


# =============================================================================
# STEP 3  Local ancestry — windowed PCA across training scenarios
# =============================================================================
# We demonstrate all three windowed PCA modes and both training scenarios:
#
#   Scenario A — IS_1KG_COHORT=FALSE (external study cohort)
#     3a. Mode A: fit per-window PCAs on reference, project study samples in
#     3b. Mode C: project new individual into saved training window PCAs
#
#   Scenario B — IS_1KG_COHORT=TRUE (1KGP3 base model)
#     3c. Mode B: fit per-window PCAs on the 1KGP3 matrix directly (no ref)
#         FLARE is invalid in this scenario — demonstrated in Step 8
#
#   3d. Unified dispatcher
cat("\n── Step 3: Windowed PCA local ancestry ──\n")

# ── 3a. IS_1KG_COHORT=FALSE — Mode A (reference-based) ──────────────────────
# Fit per-window PCAs on the reference cohort, project study samples in.
# This ensures window PC axes are stable and anchored to the reference population
# structure regardless of what study samples are present.
# Corresponds to IS_1KG_COHORT=FALSE in genoformerR_build_model.R.

local_anc_arr <- gf_windowed_pca(
  dosage_mat     = study_dosage,
  snp_chroms     = snp_chroms,
  ref_dosage_mat = ref_dosage,    # reference: PCA fitted here per window
  k_pops         = 5L,
  window_size    = 50L,           # small for toy (500 SNPs)
  step_size      = 25L
)

cat(sprintf("Local ancestry dims: [%s]\n",
            paste(dim(local_anc_arr), collapse = " x ")))

# Row-sum check — every (sample, SNP) slice must sum to 1
row_sums <- apply(local_anc_arr, c(1, 2), sum)
cat(sprintf("Row sums: min=%.6f  max=%.6f (expect ~1)\n",
            min(row_sums), max(row_sums)))

# Save window PCAs — REQUIRED for scoring new individuals later (Mode C)
wpca_rds_path <- tempfile(fileext = ".rds")
gf_save_window_pcas(local_anc_arr, wpca_rds_path)

# Summary: mean pseudo-ancestry per population per chromosome
la_summary <- gf_local_anc_summary(local_anc_arr, snp_chroms, study_pop_labels)
cat("\nLocal ancestry summary (first 10 rows):\n")
print(head(la_summary, 10))


# ── 3b. Scoring a new individual — Mode C (project into saved PCAs) ──────────
cat("\n── Mode C: scoring new individual with saved window PCAs ──\n")

# Simulate one new individual (single-sample scoring)
new_ind_dosage        <- matrix(
  vapply(pgs_weights, function(maf_proxy) {
    maf <- abs(maf_proxy) + 0.1
    p   <- min(maf, 0.5)
    sample(c(0L,1L,2L), 1L, prob=c((1-p)^2, 2*p*(1-p), p^2))
  }, integer(1L)),
  nrow = 1L
)
colnames(new_ind_dosage) <- snp_ids
rownames(new_ind_dosage) <- "new_sample_001"
storage.mode(new_ind_dosage) <- "numeric"

# Load saved window PCAs and project new individual into training subspace
ref_pcas      <- gf_load_window_pcas(wpca_rds_path)
la_new_ind    <- gf_windowed_pca(
  dosage_mat      = new_ind_dosage,
  snp_chroms      = snp_chroms,
  ref_window_pcas = ref_pcas         # Mode C: project, don't refit PCA
)
cat(sprintf("New individual local ancestry dims: [%s]\n",
            paste(dim(la_new_ind), collapse = " x ")))


# ── 3c. IS_1KG_COHORT=TRUE — Mode B (study-only, no reference) ───────────────
# When training cohort IS 1KGP3, no separate reference is needed or valid.
# PCA is fitted directly on dosage_mat (2504-sample 1KGP3; well-conditioned).
# FLARE cannot be used in this scenario — 1KGP3 cannot be both the query and
# the reference panel simultaneously (emission probability contamination).
# This corresponds to IS_1KG_COHORT=TRUE in genoformerR_build_model.R.
# Here we use ref_dosage (our 1KGP3 stand-in) as the training cohort.

cat("\n── 3c. IS_1KG_COHORT=TRUE — Mode B (1KGP3 base model, no reference) ──\n")

la_1kg_modeB <- gf_windowed_pca(
  dosage_mat  = ref_dosage,       # ref_dosage IS the training cohort here
  snp_chroms  = snp_chroms,
  # ref_dosage_mat NOT passed — dosage_mat is 1KGP3 itself (Mode B)
  k_pops      = 5L,
  window_size = 50L,
  step_size   = 25L
)
cat(sprintf("Mode B dims: [%s]\n",
            paste(dim(la_1kg_modeB), collapse = " x ")))

rs_b <- apply(la_1kg_modeB, c(1,2), sum)
cat(sprintf("Row sums: min=%.6f  max=%.6f\n", min(rs_b), max(rs_b)))

# Save Mode B window PCAs — used exactly like Mode A at scoring time (Mode C)
wpca_1kg_path <- tempfile(fileext = ".rds")
gf_save_window_pcas(la_1kg_modeB, wpca_1kg_path)
cat("Mode B window PCAs saved. Scoring new individuals uses Mode C identically.\n")

# Demonstrate Mode C on the same new individual using Mode B window PCAs
ref_pcas_1kg   <- gf_load_window_pcas(wpca_1kg_path)
la_new_modeC_1kg <- gf_windowed_pca(
  dosage_mat      = new_ind_dosage,
  snp_chroms      = snp_chroms,
  ref_window_pcas = ref_pcas_1kg  # Mode C — project into Mode B training axes
)
cat(sprintf("New individual via Mode C (from Mode B training): [%s]\n",
            paste(dim(la_new_modeC_1kg), collapse = " x ")))


# ── 3d. Unified dispatcher ───────────────────────────────────────────────────
cat("\n── 3d. Unified dispatcher (gf_local_ancestry) ──\n")
la_via_dispatcher <- gf_local_ancestry(
  method         = "windowed_pca",
  dosage_mat     = study_dosage,
  snp_chroms     = snp_chroms,
  ref_dosage_mat = ref_dosage,
  k_pops         = 5L,
  window_size    = 50L,
  step_size      = 25L
)
cat(sprintf("Dispatcher dims: [%s]\n",
            paste(dim(la_via_dispatcher), collapse = " x ")))


# =============================================================================
# STEP 4  Model build, training and inference  (requires Python / PyTorch)
# =============================================================================
cat("\n── Step 4: Python backend ──\n")

python_ok <- tryCatch({
  gf_init()
  TRUE
}, error = function(e) {
  cat("Python backend unavailable:", conditionMessage(e), "\n")
  cat("Steps 4–6 skipped. Run gf_init() after setting up the conda env.\n")
  FALSE
})

if (python_ok) {

  # ── 4a. Build model ──────────────────────────────────────────────────────────
  model <- gf_build_model(
    n_snps        = ncol(study_dosage),
    n_ld_blocks   = max(ld_blocks) + 1L,
    d_model       = 64L,
    n_heads       = 4L,
    n_layers      = 2L,
    anc_dim       = ncol(anc$anc_probs),      # 5
    k_pops        = dim(local_anc_arr)[3],    # 5
    d_local       = 16L,
    use_local_anc = TRUE,
    dropout       = 0.1
  )
  print(model)

  # Global-only comparison
  model_global <- gf_build_model(
    n_snps      = ncol(study_dosage),
    n_ld_blocks = max(ld_blocks) + 1L,
    d_model = 64L, n_heads = 4L, n_layers = 2L,
    use_local_anc = FALSE
  )
  cat(sprintf("Params: dual-level=%s  global-only=%s\n",
              format(model$n_params,        big.mark=","),
              format(model_global$n_params, big.mark=",")))

  # ── 4b. Train ────────────────────────────────────────────────────────────────
  trained <- gf_train(
    model        = model,
    dosage_mat   = study_dosage,
    anc_probs    = anc$anc_probs,
    pgs_weights  = pgs_weights,
    ld_blocks    = ld_blocks,
    chroms       = snp_chroms,
    phenotype    = phenotype,
    pop_labels   = pop_int,
    local_anc    = local_anc_arr,   # (N x M x K) from Mode A windowed PCA (Step 3a)
                                    # IS_1KG_COHORT=TRUE would pass la_1kg_modeB here
    epochs       = 10L,
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

  p_hist <- ggplot(history, aes(x = epoch)) +
    geom_line(aes(y = train_loss, colour = "Train")) +
    geom_line(aes(y = val_loss,   colour = "Val")) +
    scale_colour_manual(values = c(Train="#2E75B6", Val="#E36C09"), name=NULL) +
    labs(title="GenoFormer training loss", x="Epoch", y="Loss") +
    theme_minimal(base_size=12)
  print(p_hist)

  # ── 4c. Save and reload ──────────────────────────────────────────────────────
  model_path <- tempfile(fileext = ".pt")
  gf_save_model(model, path = model_path, save_config = TRUE)
  model_loaded <- gf_load_model(model_path)
  cat("Model saved and reloaded.\n")


  # ── Step 5: PRS inference ────────────────────────────────────────────────────
  cat("\n── Step 5: PRS computation ──\n")

  prs <- gf_predict(
    model          = model,
    dosage_mat     = study_dosage,
    anc_probs      = anc$anc_probs,
    pgs_weights    = pgs_weights,
    ld_blocks      = ld_blocks,
    chroms         = snp_chroms,
    local_anc      = local_anc_arr,
    batch_size     = 64L
  )
  prs$population <- study_pop_labels
  prs$true_pheno <- phenotype

  cat("PRS output (first 6 rows):\n")
  print(head(prs[, c("transformer_prs", "classical_prs",
                      "local_anc_available", "population")]))

  cat(sprintf("\nCorrelation with phenotype:\n"))
  cat(sprintf("  Transformer PRS: r = %.3f\n",
              cor(prs$transformer_prs, phenotype)))
  cat(sprintf("  Classical PRS:   r = %.3f\n",
              cor(prs$classical_prs,   phenotype)))


  # ── Step 5b. Score the new individual using saved artefacts ─────────────────
  cat("\n── Scoring new individual ──\n")

  # Project new individual onto training PCA (Mode A pattern)
  anc_loaded <- readRDS(anc_rds_path)
  new_pcs    <- stats::predict(anc_loaded$pca_obj,
                               newdata = new_ind_dosage)[, seq_len(ncol(anc_loaded$pcs)), drop=FALSE]
  new_anc_probs <- predict(anc_loaded$clf_obj,
                           data = as.data.frame(new_pcs),
                           type = "response")$predictions
  colnames(new_anc_probs) <- POP_NAMES

  # local ancestry already computed in Step 3b (la_new_ind)
  prs_new <- gf_predict(
    model          = model,
    dosage_mat     = new_ind_dosage,
    anc_probs      = new_anc_probs,
    pgs_weights    = pgs_weights,
    ld_blocks      = ld_blocks,
    chroms         = snp_chroms,
    local_anc      = la_new_ind,
    batch_size     = 128L
  )
  cat("New individual PRS:\n")
  print(prs_new[, c("transformer_prs", "classical_prs", "local_anc_available")])


  # ── Step 6: Portability evaluation ──────────────────────────────────────────
  cat("\n── Step 6: Portability evaluation ──\n")

  eval_res <- gf_evaluate(
    prs_tbl     = prs,
    phenotype   = phenotype,
    pop_labels  = study_pop_labels,
    ref_pop     = "EUR",
    n_bootstrap = 200L
  )
  cat("Portability metrics:\n")
  print(eval_res$metrics)
  print(eval_res$portability_plot)
  print(eval_res$distribution_plot)

}  # end python_ok


# =============================================================================
# STEP 7  Classical PRS baseline (pure R)
# =============================================================================
cat("\n── Step 7: Classical PRS baseline ──\n")

classical_prs <- scale(as.numeric(study_dosage %*% pgs_weights))[, 1]
cat("Classical PRS R² per population:\n")
for (pop in POP_NAMES) {
  mask <- study_pop_labels == pop
  cat(sprintf("  %-4s: R² = %.4f\n", pop,
              cor(classical_prs[mask], phenotype[mask])^2))
}


# =============================================================================
# STEP 8  Individual function demos
# =============================================================================
cat("\n── Step 8: Individual function demos ──\n")

# 8a. gf_windowed_pca with step_size < window_size (overlapping)
la_overlap <- gf_windowed_pca(
  dosage_mat     = study_dosage,
  snp_chroms     = snp_chroms,
  ref_dosage_mat = ref_dosage,
  k_pops         = 3L,
  window_size    = 40L,
  step_size      = 20L
)
cat(sprintf("Overlapping windows dims: [%s]\n",
            paste(dim(la_overlap), collapse=" x ")))

# 8b. gf_save_window_pcas error when attribute is missing
cat("\ngf_save_window_pcas on array without attribute (expect error):\n")
arr_no_attr <- la_overlap
attr(arr_no_attr, "window_pcas") <- NULL
tryCatch(
  gf_save_window_pcas(arr_no_attr, tempfile()),
  error = function(e) cat("Expected error:", conditionMessage(e), "\n")
)

# 8c. FLARE dispatcher rejects missing args
cat("\nFLARE missing-arg check:\n")
tryCatch(
  gf_local_ancestry(method = "flare", dosage_mat = study_dosage),
  error = function(e) cat("Expected error:", conditionMessage(e), "\n")
)

# 8d. gf_run_pipeline rejects FLARE when is_1kg_cohort=TRUE
# This demonstrates the validation guard added in v0.3.1.
# When the training cohort IS 1KGP3, FLARE cannot be used because 1KGP3
# cannot serve as its own reference panel (emission probability contamination).
cat("\ngf_run_pipeline IS_1KG_COHORT=TRUE + FLARE validation error:\n")
tryCatch(
  gf_run_pipeline(
    dosage_raw       = "data/kg3.raw",
    pgs_file         = "data/PGS000018.txt.gz",
    panel_file       = "data/kg3.panel",
    phenotype        = "data/phenotype.txt",
    is_1kg_cohort    = TRUE,
    local_anc_method = "flare"       # blocked: FLARE requires independent reference
  ),
  error = function(e) cat("Expected error:", conditionMessage(e), "\n")
)

# 8e. Internal SNP metadata helpers
cat("\nInternal SNP helpers:\n")
test_ids <- c("chr1_1000_A_T", "chr3_5000_G_C", "chr5_99000_T_A")
cat("  chromosomes:", paste(genoformerR:::.snp_chroms(test_ids),    collapse=" "), "\n")
cat("  positions:  ", paste(genoformerR:::.snp_positions(test_ids), collapse=" "), "\n")
cat("  LD blocks:  ", paste(genoformerR:::.assign_ld_blocks(
  test_ids, genoformerR:::.snp_chroms(test_ids)), collapse=" "), "\n")


# =============================================================================
# SESSION INFO
# =============================================================================
cat("\n── Session info ──\n")
cat(sprintf("R: %s\n", R.version$version.string))
cat(sprintf("genoformerR: %s\n", as.character(packageVersion("genoformerR"))))
if (exists("python_ok") && python_ok) {
  torch_ver <- tryCatch(
    reticulate::py_to_r(reticulate::import("torch")$`__version__`),
    error = function(e) "unavailable"
  )
  cat(sprintf("PyTorch: %s\n", torch_ver))
}
cat("\nToy example complete.\n")
