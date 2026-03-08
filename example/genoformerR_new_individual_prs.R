# =============================================================================
# genoformerR — Ancestry-Aware PRS for a New Individual (v0.3.1)
# =============================================================================
#
# Scores a new VCF file using a previously trained GenoFormer model.
# Handles both global and local ancestry inference consistently with training.
#
# HOW TRAINING SCENARIO AFFECTS SCORING
# ─────────────────────────────────────────────────────────────────────────────
# The training script (genoformerR_build_model.R) has an IS_1KG_COHORT flag
# that controls how ancestry was inferred during training. This determines what
# methods are valid here at scoring time.
#
# IS_1KG_COHORT = TRUE (1KGP3 base model)
#   Training used windowed PCA Mode B (no external reference).
#   LOCAL_ANC_METHOD for scoring: "windowed_pca" or "flare".
#     - "windowed_pca": uses Mode C — project into window_pcas.rds. The
#       LocalAncestryProjector MLP was trained on windowed PCA K-vectors,
#       so scoring must produce the same kind of input.
#     - "flare": valid at scoring time (new individual is never in 1KGP3),
#       BUT only if training also used FLARE. If training used windowed_pca,
#       you must score with windowed_pca — the MLP learned from windowed PCA
#       pseudo-ancestry dimensions, not FLARE's named-population posteriors.
#
# IS_1KG_COHORT = FALSE (external study cohort)
#   Training used windowed PCA Mode A (projected into 1KGP3) or FLARE.
#   LOCAL_ANC_METHOD for scoring must match training exactly.
#     - "flare" training → "flare" scoring with the same ref_vcf + ref_panel
#     - "windowed_pca" training → "windowed_pca" scoring with window_pcas.rds
#
# In all cases: set LOCAL_ANC_METHOD below to match what was used in training.
#
# PREREQUISITES — files produced during training (keep alongside the model)
# ─────────────────────────────────────────────────────────────────────────────
#   genoformer_best.pt            trained model weights
#   genoformer_best_config.json   architecture config (written by gf_save_model)
#   ancestry.rds                  gf_ancestry object — contains pca_obj + clf_obj
#   window_pcas.rds               per-window prcomp objects (windowed PCA only)
#   pgs_snps.txt                  SNP IDs used during training
#   PGS000018.txt.gz              PGS Catalog scoring file
#   kg3_ref.vcf.gz                1KGP3 reference VCF (FLARE method only)
#   ref_panel.txt                 sample_id <TAB> population (FLARE only)
#   genetic_maps/                 per-chromosome recombination maps (FLARE only)
#   train_prs_scores.rds          training PRS output (for percentile reporting)
#
# NEW INDIVIDUAL FILE
# ─────────────────────────────────────────────────────────────────────────────
#   new_individual.vcf.gz         single-sample or multi-sample VCF (unphased OK)
#
# =============================================================================

library(genoformerR)
library(data.table)

# ── Paths ─────────────────────────────────────────────────────────────────────
MODEL_PT     <- "models/genoformer_best.pt"
MODEL_CFG    <- "models/genoformer_best_config.json"
ANCESTRY_RDS <- "models/ancestry.rds"
WINDOW_PCAS  <- "models/window_pcas.rds"       # only needed for windowed_pca method
PGS_FILE     <- "data/PGS000018.txt.gz"
PGS_SNPS     <- "data/pgs_snps.txt"
NEW_VCF      <- "new_individual.vcf.gz"
REF_VCF      <- "data/kg3_ref.vcf.gz"          # FLARE only
REF_PANEL    <- "data/ref_panel.txt"            # FLARE only: sample_id <TAB> pop
GMAP_DIR     <- "data/genetic_maps"             # FLARE only
FLARE_JAR    <- "tools/flare.jar"               # FLARE only
OUTDIR       <- "results/new_individual"
PLINK2       <- "plink2"

# LOCAL ANCESTRY METHOD — must match what was used during training.
#
#   "flare"        Valid when training used FLARE. Requires Java + 1KGP3 ref VCF.
#                  Also valid as a first choice when IS_1KG_COHORT=TRUE was used
#                  during training with windowed_pca — but only if training itself
#                  used FLARE. Never mix methods between training and scoring.
#
#   "windowed_pca" Required when training used windowed_pca (IS_1KG_COHORT=TRUE
#                  Mode B, or IS_1KG_COHORT=FALSE Mode A). Uses saved window_pcas.rds
#                  to project into the training PC subspace (Mode C).
LOCAL_ANC_METHOD <- "flare"   # set to match training: "flare" | "windowed_pca"

fs::dir_create(OUTDIR)


# =============================================================================
# STEP 1  Initialise Python backend
# =============================================================================
gf_init()


# =============================================================================
# STEP 2  Extract PGS SNPs from the new VCF
# =============================================================================
# Restrict to exactly the SNP set used during training.

system2(PLINK2, c(
  "--vcf",         NEW_VCF,
  "--make-pgen",
  "--max-alleles", "2",
  "--rm-dup",      "exclude-all",
  "--out",         file.path(OUTDIR, "new_ind_raw"),
  "--silent"
))

system2(PLINK2, c(
  "--pfile",  file.path(OUTDIR, "new_ind_raw"),
  "--extract", PGS_SNPS,
  "--export", "A",
  "--out",    file.path(OUTDIR, "new_ind_dosage"),
  "--silent"
))

raw        <- data.table::fread(
  file.path(OUTDIR, "new_ind_dosage.raw"), nThread = 4L
)
sample_ids <- raw$IID
dosage_mat <- as.matrix(raw[, -(1:6)])
dosage_mat[is.na(dosage_mat)] <- 0
storage.mode(dosage_mat) <- "numeric"

cat(sprintf("Dosage matrix: %d sample(s) x %d SNPs\n",
            nrow(dosage_mat), ncol(dosage_mat)))


# =============================================================================
# STEP 3  Reconstruct SNP metadata
# =============================================================================
snp_ids       <- colnames(dosage_mat)
snp_chroms    <- genoformerR:::.snp_chroms(snp_ids)
snp_positions <- genoformerR:::.snp_positions(snp_ids)
ld_blocks     <- genoformerR:::.assign_ld_blocks(snp_ids, snp_chroms)

# If SNP IDs are plain rsIDs (no embedded chr:pos), use the training .pvar:
#   pvar          <- data.table::fread("data/plink/kg3_qc.pvar")
#   pos_map       <- pvar[match(snp_ids, pvar$ID)]
#   snp_chroms    <- as.integer(sub("chr","", pos_map$`#CHROM`))
#   snp_positions <- pos_map$POS
#   ld_blocks     <- genoformerR:::.assign_ld_blocks(snp_ids, snp_chroms)


# =============================================================================
# STEP 4  Load PGS effect weights
# =============================================================================
pgs        <- data.table::fread(PGS_FILE, skip = "#", nThread = 2L)
names(pgs) <- tolower(names(pgs))
id_col     <- intersect(c("rsid","variant_id","id","snp"),        names(pgs))[1]
w_col      <- intersect(c("effect_weight","beta","or"),            names(pgs))[1]
ea_col     <- intersect(c("effect_allele","a1"),                   names(pgs))[1]

pgs_matched    <- pgs[match(snp_ids, pgs[[id_col]]), ]
pgs_weights    <- as.numeric(pgs_matched[[w_col]])
effect_alleles <- as.character(pgs_matched[[ea_col]])

n_miss <- sum(is.na(pgs_weights))
if (n_miss > 0)
  warning(sprintf("%d SNPs had no matching weight — set to 0", n_miss))
pgs_weights[is.na(pgs_weights)] <- 0


# =============================================================================
# STEP 5  Global ancestry — project into training reference space
# =============================================================================
# Load the saved gf_ancestry object from training.
# The pca_obj inside was fitted either on 1KGP3 directly (Mode B, when training
# used IS_1KG_COHORT=TRUE) or on 1KGP3 as a reference with study samples
# projected in (Mode A, IS_1KG_COHORT=FALSE). In both cases it is a prcomp
# object whose rotation matrix defines a stable PC subspace.
#
# stats::predict.prcomp applies the same centering vector and rotation to the
# new individual's dosage values, placing them in the identical PC space used
# during training. This works correctly for both Mode A and Mode B pca_obj.
#
# Never refit PCA on the new individual alone — a 1-sample PCA has no variance
# to decompose and produces axes unrelated to the training subspace.

anc_train <- readRDS(ANCESTRY_RDS)

# Project new individual(s) into training PC space
new_pcs <- stats::predict(
  anc_train$pca_obj,
  newdata = dosage_mat
)[, seq_len(ncol(anc_train$pcs)), drop = FALSE]

# Predict ancestry probabilities using the trained Random Forest
if (inherits(anc_train$clf_obj, "ranger")) {
  anc_probs <- predict(
    anc_train$clf_obj,
    data = as.data.frame(new_pcs),
    type = "response"
  )$predictions
} else {
  # Python sklearn GBM backend
  anc_probs <- reticulate::py_to_r(
    anc_train$clf_obj$predict_proba(reticulate::r_to_py(new_pcs))
  )
}
colnames(anc_probs) <- c("AFR","AMR","EAS","EUR","SAS")

inferred_pop <- colnames(anc_probs)[apply(anc_probs, 1, which.max)]
entropy      <- -rowSums(anc_probs * log(anc_probs + 1e-9))

cat("\nGlobal ancestry probabilities:\n")
print(round(anc_probs, 3))
cat(sprintf("Inferred population: %s (entropy = %.3f)\n",
            inferred_pop[1], entropy[1]))


# =============================================================================
# STEP 6  Local ancestry
# =============================================================================
# Choose method to match the method used during training.
# The model's LocalAncestryProjector MLP was trained on a specific kind of
# local_anc input; mixing methods between training and scoring degrades accuracy.

if (LOCAL_ANC_METHOD == "flare") {

  # ── FLARE ───────────────────────────────────────────────────────────────────
  # A new individual is never in 1KGP3, so using 1KGP3 as the FLARE reference
  # is always valid here — there is no emission-probability contamination.
  #
  # When to use FLARE here:
  #   - Training used FLARE (IS_1KG_COHORT=FALSE + local_anc_method="flare")
  #     → use same ref_vcf and ref_panel; FLARE's K population axes are
  #       anchored to reference populations, consistency is automatic.
  #   - Training used IS_1KG_COHORT=TRUE + local_anc_method="windowed_pca"
  #     → do NOT use FLARE here; the model's LocalAncestryProjector MLP was
  #       trained on windowed PCA K-vectors, not FLARE's named-population
  #       posteriors. Set LOCAL_ANC_METHOD="windowed_pca" instead.
  #
  # The same ref_vcf and ref_panel used during training must be supplied so
  # that K and population ordering are identical between training and scoring.

  flare_out <- gf_local_ancestry_flare(
    vcf_file   = NEW_VCF,
    ref_vcf    = REF_VCF,
    sample_map = REF_PANEL,
    gmap_dir   = GMAP_DIR,
    outdir     = file.path(OUTDIR, "flare"),
    chroms     = 1:22,
    flare_jar  = FLARE_JAR,
    nthreads   = 4L,
    em_its     = 3L
  )

  local_anc_arr <- gf_parse_flare(
    flare_prefixes    = flare_out,
    snp_ids           = snp_ids,
    snp_positions     = snp_positions,
    snp_chroms        = snp_chroms,
    sample_ids        = sample_ids,
    global_anc        = anc_probs,
    entropy_threshold = 0.3
  )

} else if (LOCAL_ANC_METHOD == "windowed_pca") {

  # ── Windowed PCA — Mode C (project into saved training window PCAs) ─────────
  # Load the per-window prcomp objects saved during training.
  # gf_windowed_pca projects the new individual into the same PC subspace used
  # at training time via stats::predict.prcomp per window (Mode C).
  #
  # This is correct regardless of which training mode produced the window PCAs:
  #   - IS_1KG_COHORT=TRUE (Mode B): PCAs were fitted on the 1KGP3 matrix.
  #     Mode C reprojects into those 1KGP3-derived axes.
  #   - IS_1KG_COHORT=FALSE (Mode A): PCAs were fitted on 1KGP3 as reference,
  #     with study samples projected in. Mode C reprojects into those same axes.
  # In both cases, window_pcas.rds is the bridge between training and scoring.
  #
  # Never call gf_windowed_pca without ref_window_pcas= at scoring time.
  # Refitting per-window PCAs on a single sample produces degenerate output
  # (zero variance per window, arbitrary axis orientations).

  if (!file.exists(WINDOW_PCAS))
    stop(paste(
      "window_pcas.rds not found at:", WINDOW_PCAS,
      "\nThis file must be produced during training by calling:",
      "gf_save_window_pcas(local_anc_arr, 'models/window_pcas.rds')"
    ))

  ref_pcas <- gf_load_window_pcas(WINDOW_PCAS)

  local_anc_arr <- gf_windowed_pca(
    dosage_mat      = dosage_mat,
    snp_chroms      = snp_chroms,
    ref_window_pcas = ref_pcas      # Mode C: project, do not refit PCA
  )

} else {
  stop(sprintf("Unknown LOCAL_ANC_METHOD: '%s'. Use 'flare' or 'windowed_pca'.",
               LOCAL_ANC_METHOD))
}

cat(sprintf("\nLocal ancestry array: [%s]\n",
            paste(dim(local_anc_arr), collapse = " x ")))

# Per-chromosome summary (first 3 chromosomes, first sample)
for (chr in 1:min(3L, max(snp_chroms))) {
  mask     <- snp_chroms == chr
  if (!any(mask)) next
  mean_anc <- colMeans(local_anc_arr[1, mask, , drop=FALSE])
  cat(sprintf("chr%d mean local ancestry: %s\n", chr,
              paste(names(mean_anc), round(mean_anc, 3), sep="=", collapse="  ")))
}


# =============================================================================
# STEP 7  Load the trained GenoFormer model
# =============================================================================
model <- gf_load_model(
  path   = MODEL_PT,
  config = MODEL_CFG,
  device = "auto"
)
print(model)

# Verify K alignment
if (dim(local_anc_arr)[3] != model$config$k_pops)
  stop(sprintf(
    "local_anc K=%d does not match model k_pops=%d.\n%s",
    dim(local_anc_arr)[3], model$config$k_pops,
    "Ensure the same local ancestry method and reference panel were used at training."
  ))


# =============================================================================
# STEP 8  Compute PRS
# =============================================================================
prs_result <- gf_predict(
  model          = model,
  dosage_mat     = dosage_mat,
  anc_probs      = anc_probs,
  pgs_weights    = pgs_weights,
  ld_blocks      = ld_blocks,
  chroms         = snp_chroms,
  effect_alleles = effect_alleles,
  local_anc      = local_anc_arr,
  batch_size     = 128L
)


# =============================================================================
# STEP 9  Interpret and report
# =============================================================================
prs_result$sample_id         <- sample_ids
prs_result$inferred_pop      <- inferred_pop
prs_result$admixture_entropy <- entropy

# Population-normalised percentile using training distribution as reference
# Requires train_prs_scores.rds saved at the end of the training run
if (file.exists("models/train_prs_scores.rds")) {
  train_prs  <- readRDS("models/train_prs_scores.rds")
  pop_stats  <- aggregate(transformer_prs ~ population, data = train_prs,
                           FUN = function(x) c(mean=mean(x), sd=sd(x)))

  get_percentile <- function(raw_z, pop) {
    row <- pop_stats[pop_stats$population == pop, ]
    if (nrow(row) == 0) return(NA_real_)
    pnorm((raw_z - row$transformer_prs["mean"]) /
           row$transformer_prs["sd"]) * 100
  }

  prs_result$percentile <- mapply(
    get_percentile,
    prs_result$transformer_prs,
    prs_result$inferred_pop
  )
} else {
  prs_result$percentile <- NA_real_
  message("train_prs_scores.rds not found — percentile column set to NA.")
}

cat("\n══════════════════════════════════════════════\n")
cat("  GenoFormer PRS Report\n")
cat("══════════════════════════════════════════════\n")
for (i in seq_len(nrow(prs_result))) {
  cat(sprintf("  Sample:                %s\n",  prs_result$sample_id[i]))
  cat(sprintf("  Local anc method:      %s\n",  LOCAL_ANC_METHOD))
  cat(sprintf("  Inferred population:   %s\n",  prs_result$inferred_pop[i]))
  cat(sprintf("  Admixture entropy:     %.3f\n", prs_result$admixture_entropy[i]))
  cat(sprintf("  Transformer PRS (z):   %.4f\n", prs_result$transformer_prs[i]))
  cat(sprintf("  Classical PRS (z):     %.4f\n", prs_result$classical_prs[i]))
  cat(sprintf("  Local anc available:   %s\n",  prs_result$local_anc_available[i]))
  if (!is.na(prs_result$percentile[i]))
    cat(sprintf("  Pop percentile:        %.1f%%\n", prs_result$percentile[i]))
  cat("──────────────────────────────────────────────\n")
}

write.csv(prs_result,
          file.path(OUTDIR, "prs_scores_new_individual.csv"),
          row.names = FALSE)
cat(sprintf("Saved to: %s\n",
            file.path(OUTDIR, "prs_scores_new_individual.csv")))


# =============================================================================
# NOTES ON OUTPUTS
# =============================================================================
#
# transformer_prs
#   Primary output. Z-scored within the current batch. Conditioned on both
#   global ancestry (FiLM layers) and local ancestry (token fusion), making it
#   more portable across populations than the classical score.
#
# classical_prs
#   Weighted dosage sum (C+T). Included for benchmarking. Ignores ancestry
#   context; typically shows reduced portability for admixed individuals.
#
# admixture_entropy
#   Shannon entropy of soft ancestry probabilities. Near 0 = homogeneous
#   continental ancestry; above ~0.5 = substantial admixture. Local ancestry
#   conditioning is most important for high-entropy (admixed) individuals.
#
# local_anc_available
#   FALSE if local_anc was NULL and the model fell back to broadcasting global
#   ancestry per-SNP. For admixed individuals always ensure this is TRUE.
#
#
# CONSISTENCY REQUIREMENTS — training vs scoring
# ─────────────────────────────────────────────────────────────────────────────
#
# Global ancestry
#   Always use the saved pca_obj and clf_obj from ancestry.rds.
#   stats::predict.prcomp applies the training-time centering and rotation to
#   new dosage values. This is correct whether the training pca_obj came from
#   Mode A (1KGP3 reference projection) or Mode B (1KGP3 study-only PCA) —
#   both produce a valid prcomp rotation matrix. Never refit PCA on a new
#   individual alone.
#
# Local ancestry — windowed PCA
#   Always use ref_window_pcas= with the saved window_pcas.rds (Mode C).
#   Never call gf_windowed_pca() without this argument at scoring time.
#   This applies whether training used Mode A (IS_1KG_COHORT=FALSE) or Mode B
#   (IS_1KG_COHORT=TRUE) — the window_pcas.rds bridges both to Mode C.
#
# Local ancestry — FLARE
#   Use the same ref_vcf and ref_panel that were used during training.
#   FLARE's K dimensions are anchored to named reference populations, so
#   consistency is automatic as long as the reference is unchanged.
#   FLARE cannot be used at scoring time if training used windowed_pca: the
#   LocalAncestryProjector MLP was trained on windowed PCA K-vectors and
#   expects inputs of that form, not FLARE's named-population posteriors.
#
# Method matching
#   LOCAL_ANC_METHOD here must be identical to the method used during training.
#   Mixing methods produces inputs the MLP was not trained on and silently
#   degrades the transformer's ancestry conditioning.
#   Training scenario → scoring method:
#     IS_1KG_COHORT=TRUE  + windowed_pca → windowed_pca (Mode C)
#     IS_1KG_COHORT=FALSE + windowed_pca → windowed_pca (Mode C)
#     IS_1KG_COHORT=FALSE + flare        → flare (same ref_vcf + ref_panel)
#     IS_1KG_COHORT=FALSE + rfmix2       → rfmix2 (same ref_vcf + ref_panel)
#     IS_1KG_COHORT=TRUE  + none         → none (local_anc=NULL)
#     IS_1KG_COHORT=FALSE + none         → none (local_anc=NULL)
