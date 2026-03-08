# =============================================================================
# genoformerR — Model Generation Script (v0.3.1)
# =============================================================================
#
# Full pipeline from raw QC'd genotypes to a saved, evaluated GenoFormer model.
# Produces all artefacts needed for scoring new individuals.
#
# COHORT FLAG (IS_1KG_COHORT) — READ THIS FIRST
# ─────────────────────────────────────────────────────────────────────────────
# This script supports two training scenarios which differ in how Phases 5
# and 6 behave. Set IS_1KG_COHORT in the CONFIGURATION block below.
#
#   IS_1KG_COHORT = TRUE   The genotype matrix being trained on IS the 1KGP3
#                          cohort itself (2504 samples, 26 populations).
#                          This is the standard base-model scenario: you train
#                          on 1KGP3 with simulated or real phenotypes, then
#                          apply the model to external cohorts at scoring time.
#
#                          Consequences:
#                          - Global ancestry (Phase 5): Mode B. Sample IDs in
#                            dosage_mat match the panel file directly, so PCA
#                            is fitted and the RF is trained on the same
#                            labelled samples. No separate reference needed.
#                          - Local ancestry (Phase 6): FLARE is NOT available
#                            because 1KGP3 cannot serve as both the query and
#                            the reference panel simultaneously (see note below).
#                            Use "windowed_pca" (Mode B) or "none".
#
#   IS_1KG_COHORT = FALSE  The genotype matrix is an external study cohort
#                          (e.g. a GWAS cohort, biobank, or disease study) that
#                          is NOT part of 1KGP3. 1KGP3 acts as the reference.
#
#                          Consequences:
#                          - Global ancestry (Phase 5): Mode A. PCA is fitted
#                            on 1KGP3 and study samples are projected in.
#                          - Local ancestry (Phase 6): All three methods
#                            available. FLARE uses 1KGP3 as a proper
#                            independent reference. Windowed PCA uses Mode A
#                            (1KGP3 fits per-window PCAs, study projected in).
#
# WHY FLARE CANNOT BE USED WHEN IS_1KG_COHORT = TRUE
# ─────────────────────────────────────────────────────────────────────────────
# FLARE is a Hidden Markov Model where:
#   - Emission probabilities = per-population allele frequencies from ref_vcf
#   - Transition probabilities = recombination rates from the genetic map
#
# When study and reference are the same individuals, FLARE's allele frequency
# estimates are derived from the very samples being classified. The query
# alleles contribute to the reference frequencies, biasing every emission
# probability. The resulting local ancestry posteriors are not valid.
#
# The correct alternative for 1KGP3 training is windowed PCA Mode B:
# PCA is fit per-window on the 2504 1KGP3 samples (large, population-diverse,
# well-conditioned), producing K pseudo-ancestry dimensions that capture
# within-chromosome population structure without needing an external reference.
# The per-window prcomp objects are saved to window_pcas.rds and used at
# scoring time to project new individuals into the same subspace.
#
# INPUT FILES REQUIRED
# ─────────────────────────────────────────────────────────────────────────────
#   study.vcf.gz             study/training cohort VCF (unphased OK)
#   PGS000018.txt.gz         PGS Catalog scoring file (.txt.gz)
#   phenotype.txt            2-col TSV: sample_id <TAB> phenotype_value
#
#   1KGP3 reference (download via gf_download_1kg or provide existing):
#     data/1kg/ALL.chr*.vcf.gz
#     data/1kg/*.panel
#
#   FLARE only (IS_1KG_COHORT = FALSE and LOCAL_ANC_METHOD = "flare"):
#     data/ref_panel.txt             2-col: sample_id <TAB> population
#     data/genetic_maps/             per-chromosome recombination maps
#     tools/flare.jar                Java >= 11 required
#
# OUTPUT FILES (all in OUTDIR)
# ─────────────────────────────────────────────────────────────────────────────
#   genoformer_best.pt          model weights
#   genoformer_best_config.json architecture config
#   ancestry.rds                PCA + RF objects for global ancestry projection
#   window_pcas.rds             per-window prcomp objects (windowed_pca only)
#   pgs_snps.txt                SNP list for new-individual PLINK2 extraction
#   train_prs_scores.rds        training PRS + pop labels for percentile norms
#   training_history.csv        loss curve
#   portability.png             R² portability plot across superpopulations
#
# LOCAL ANCESTRY METHODS — compatibility matrix
# ─────────────────────────────────────────────────────────────────────────────
#   Method          IS_1KG_COHORT=TRUE   IS_1KG_COHORT=FALSE   Notes
#   "flare"         NOT VALID            Valid                  Needs Java + ref VCF
#   "windowed_pca"  Valid (Mode B)       Valid (Mode A)         No external tools
#   "none"          Valid                Valid                  Global-only baseline
# =============================================================================

library(genoformerR)
library(data.table)
library(ggplot2)


# =============================================================================
# CONFIGURATION — edit this block
# =============================================================================

# -- Input files ---------------------------------------------------------------
STUDY_VCF    <- "data/study.vcf.gz"          # study cohort VCF
PGS_FILE     <- "data/PGS000018.txt.gz"      # PGS Catalog file
PHENOTYPE    <- "data/phenotype.txt"          # 2-col: sample_id  phenotype
PANEL_FILE   <- "data/1kg/integrated_call_samples_v3.20130502.ALL.panel"
KG3_VCF_DIR  <- "data/1kg"                   # directory of 1KGP3 per-chr VCFs

# -- FLARE inputs (only used when LOCAL_ANC_METHOD = "flare") ------------------
REF_PANEL    <- "data/ref_panel.txt"          # sample_id <TAB> population
GMAP_DIR     <- "data/genetic_maps"           # recombination maps
FLARE_JAR    <- "tools/flare.jar"

# -- Outputs -------------------------------------------------------------------
OUTDIR       <- "models/genoformer_v1"
PLINK2       <- "plink2"

# -- Method --------------------------------------------------------------------
# IS_1KG_COHORT: TRUE  = training cohort IS 1KGP3 (base model scenario)
#                FALSE = training cohort is an external study cohort
IS_1KG_COHORT    <- TRUE

LOCAL_ANC_METHOD <- "windowed_pca"  # "windowed_pca" | "none"
                                    # "flare" requires IS_1KG_COHORT = FALSE
CHROMS           <- 1:22            # set to 21:22 for a quick test run

# -- Architecture (override if needed) ----------------------------------------
D_MODEL      <- 256L
N_HEADS      <- 8L
N_LAYERS     <- 6L
D_LOCAL      <- 64L
K_POPS       <- 5L    # 5 superpopulations; match ref_panel populations

# -- Training ------------------------------------------------------------------
EPOCHS       <- 100L
BATCH_SIZE   <- 64L
LR           <- 1e-4
WEIGHT_DECAY <- 0.01
ANC_WEIGHT   <- 0.3   # weight on population CE auxiliary loss
CALIB_WEIGHT <- 0.5   # weight on ancestry calibration loss
VAL_SPLIT    <- 0.1

set.seed(42L)
fs::dir_create(OUTDIR)


# =============================================================================
# PHASE 0  Setup
# =============================================================================

cat("════════════════════════════════════════════════\n")
cat(" GenoFormer Model Generation\n")
cat(sprintf(" Training cohort:  %s\n",
            if (IS_1KG_COHORT) "1KGP3 (base model)" else "External study cohort"))
cat(sprintf(" Local ancestry:   %s\n", LOCAL_ANC_METHOD))
cat(sprintf(" Output:           %s\n", OUTDIR))
cat("════════════════════════════════════════════════\n\n")

# Validate method / cohort combination
if (IS_1KG_COHORT && LOCAL_ANC_METHOD == "flare")
  stop(paste(
    "LOCAL_ANC_METHOD='flare' is incompatible with IS_1KG_COHORT=TRUE.",
    "\n  FLARE requires an independent reference panel.",
    "\n  When the training cohort IS 1KGP3, using 1KGP3 as both the query",
    "\n  and the reference contaminates allele frequency estimates in FLARE's",
    "\n  emission probabilities, producing invalid local ancestry posteriors.",
    "\n  Use LOCAL_ANC_METHOD='windowed_pca' or 'none' for 1KGP3 base model training."
  ))

# Initialise Python / PyTorch backend
# Creates conda env "genoformer" on first run (~5 min); instant on subsequent runs
gf_init()


# =============================================================================
# PHASE 1  Download 1KGP3 reference (skip if already downloaded)
# =============================================================================
cat("[1/8] Checking 1KGP3 reference data...\n")

kg3_files <- gf_download_1kg(
  outdir     = KG3_VCF_DIR,
  chroms     = CHROMS,
  panel_only = FALSE,       # set TRUE to skip VCF download (windowed_pca + none)
  ncores     = 4L
)
# After download:
#   data/1kg/ALL.chr{N}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
#   data/1kg/integrated_call_samples_v3.20130502.ALL.panel


# =============================================================================
# PHASE 2  Genotype QC and dosage extraction
# =============================================================================
cat("[2/8] Genotype QC...\n")

qc <- gf_qc(
  vcf_dir   = KG3_VCF_DIR,    # directory with per-chromosome VCFs
  pgs_file  = PGS_FILE,
  outdir    = file.path(OUTDIR, "plink"),
  maf       = 0.01,
  hwe       = 1e-6,
  geno_miss = 0.02,
  mind_miss = 0.02,
  ld_window = 200L,
  ld_step   = 50L,
  ld_r2     = 0.3,
  plink2    = PLINK2,
  chroms    = CHROMS
)
cat(sprintf("QC complete: %d samples  %d PGS SNPs\n", qc$n_samples, qc$n_snps_pgs))

# Save the PGS SNP list — needed to extract matching SNPs from new individual VCFs
file.copy("data/plink/pgs_snps.txt",
          file.path(OUTDIR, "pgs_snps.txt"), overwrite = TRUE)

# Load the dosage matrix
raw        <- data.table::fread(qc$dosage_raw, nThread = 4L)
sample_ids <- raw$IID
dosage_mat <- as.matrix(raw[, -(1:6)])
dosage_mat[is.na(dosage_mat)] <- 0
storage.mode(dosage_mat) <- "numeric"

N <- nrow(dosage_mat)
M <- ncol(dosage_mat)
cat(sprintf("Dosage matrix: %d x %d\n", N, M))

# Derive SNP metadata from column names
snp_ids       <- colnames(dosage_mat)
snp_chroms    <- genoformerR:::.snp_chroms(snp_ids)
snp_positions <- genoformerR:::.snp_positions(snp_ids)
ld_blocks     <- genoformerR:::.assign_ld_blocks(snp_ids, snp_chroms)
cat(sprintf("LD blocks: %d\n", max(ld_blocks) + 1L))


# =============================================================================
# PHASE 3  PGS effect weights
# =============================================================================
cat("[3/8] Loading PGS weights...\n")

pgs        <- data.table::fread(PGS_FILE, skip = "#", nThread = 2L)
names(pgs) <- tolower(names(pgs))
id_col     <- intersect(c("rsid","variant_id","id","snp"),  names(pgs))[1]
w_col      <- intersect(c("effect_weight","beta","or"),     names(pgs))[1]
ea_col     <- intersect(c("effect_allele","a1"),            names(pgs))[1]

pgs_matched    <- pgs[match(snp_ids, pgs[[id_col]]), ]
pgs_weights    <- as.numeric(pgs_matched[[w_col]])
effect_alleles <- as.character(pgs_matched[[ea_col]])
pgs_weights[is.na(pgs_weights)] <- 0

n_zero <- sum(pgs_weights == 0)
cat(sprintf("PGS weights loaded: %d non-zero  %d zero-weight\n",
            M - n_zero, n_zero))


# =============================================================================
# PHASE 4  Phenotype
# =============================================================================
cat("[4/8] Loading phenotype...\n")

ph_df     <- data.table::fread(PHENOTYPE, nThread = 2L)
phenotype <- as.numeric(ph_df[[2]][match(sample_ids, ph_df[[1]])])
n_missing <- sum(is.na(phenotype))
if (n_missing > 0)
  warning(sprintf("%d samples have no phenotype — they will be excluded from training", n_missing))
cat(sprintf("Phenotype: N=%d  mean=%.3f  sd=%.3f  missing=%d\n",
            sum(!is.na(phenotype)), mean(phenotype, na.rm=TRUE),
            sd(phenotype, na.rm=TRUE), n_missing))


# =============================================================================
# PHASE 5  Global ancestry inference
# =============================================================================
# IS_1KG_COHORT = TRUE  → Mode B: dosage_mat IS 1KGP3, so sample IDs match
#                          the panel file directly. PCA is fitted on dosage_mat
#                          and the RF classifier is trained on the same samples'
#                          known labels. No separate reference needed or valid.
#
# IS_1KG_COHORT = FALSE → Mode A: dosage_mat is an external cohort whose IDs
#                          do not appear in the panel file. PCA is fitted on the
#                          1KGP3 QC dosage (ref_dosage_raw=) and external samples
#                          are projected into that fixed PC subspace. The RF is
#                          trained on 1KGP3 labels and predicts for external samples.
cat("[5/8] Global ancestry inference...\n")

anc <- gf_ancestry(
  dosage_mat     = dosage_mat,
  panel_file     = PANEL_FILE,
  ref_dosage_raw = if (IS_1KG_COHORT) NULL else qc$dosage_raw,
  n_pcs          = 20L,
  method         = "r",
  n_trees        = 200L,
  seed           = 42L
)
cat(sprintf("Mode %s ancestry inference.\n",
            if (IS_1KG_COHORT) "B (study-only)" else "A (reference projection)"))
print(anc)

# Save the full object — pca_obj and clf_obj inside are essential for scoring
saveRDS(anc, file.path(OUTDIR, "ancestry.rds"))
cat(sprintf("ancestry.rds saved.\n"))

# PCA plot
p_pca <- plot(anc)
ggplot2::ggsave(file.path(OUTDIR, "ancestry_pca.png"), p_pca,
                width = 8, height = 6, dpi = 150)

# Integer-coded population labels for the training CE loss
pop_int <- as.integer(factor(anc$labels,
                              levels = c("AFR","AMR","EAS","EUR","SAS"))) - 1L


# =============================================================================
# PHASE 6  Local ancestry
# =============================================================================
# The local_anc array produced here is passed directly to gf_train() as the
# local_anc= argument. Its shape is (N x M x K): one K-vector of ancestry
# posteriors per sample per SNP. The model's LocalAncestryProjector MLP reads
# these K-vectors at every SNP position and fuses them into the SNP token
# embeddings before attention.
#
# IS_1KG_COHORT = TRUE — windowed PCA Mode B
# ─────────────────────────────────────────────────────────────────────────────
# FLARE is unavailable: 1KGP3 cannot be its own FLARE reference panel.
# Windowed PCA does not need a reference panel; it fits PCA per genomic window
# directly on the 2504-sample 1KGP3 dosage matrix (N is large, population-
# diverse, per-window PCA is well-conditioned). Output is (N x M x K_POPS)
# where K = k_pops pseudo-ancestry dimensions (softmax-normalised PC scores).
#
# IMPORTANT: gf_save_window_pcas() saves the per-window prcomp objects.
# These MUST be preserved. At scoring time, new individuals are projected
# into the training window PC subspace using gf_load_window_pcas() with
# ref_window_pcas= (Mode C). Without this, each new individual's window PCA
# would be fit from scratch — degenerate for a single sample (all variances
# are zero) and incomparable to the training subspace even with more samples.
#
# IS_1KG_COHORT = FALSE — method-dependent
# ─────────────────────────────────────────────────────────────────────────────
# FLARE: Uses 1KGP3 as an independent reference panel. Valid because study
#   samples are absent from the reference, so allele frequency estimates are
#   uncontaminated. FLARE's K axes are named biological populations, so the
#   scoring coordinate system is inherently consistent — no equivalent of
#   gf_save_window_pcas is needed.
#
# Windowed PCA (Mode A): PCA fitted on 1KGP3 per window, study samples
#   projected in via ref_dosage_mat=. Window PCAs must still be saved.
#
# Local ancestry = NULL (method "none"):
# ─────────────────────────────────────────────────────────────────────────────
# local_anc_arr remains NULL. gf_build_model() is called with use_local_anc=FALSE
# and gf_train() receives local_anc=NULL. The model trains using only global
# ancestry conditioning (FiLM layers). No LocalAncestryProjector MLP is built.
cat(sprintf("[6/8] Local ancestry — %s...\n", LOCAL_ANC_METHOD))

local_anc_arr <- NULL  # stays NULL for LOCAL_ANC_METHOD = "none"

if (LOCAL_ANC_METHOD == "windowed_pca") {

  if (IS_1KG_COHORT) {

    # ── 1KGP3 base model — Mode B ──────────────────────────────────────────────
    # dosage_mat IS the 1KGP3 cohort. PCA is fitted on dosage_mat per window
    # with no separate reference. The 2504 samples spanning 5 superpopulations
    # provide sufficient variance for stable per-window PCs.
    # Do NOT pass ref_dosage_mat= here: that would load the same file twice
    # (dosage_mat and ref_dosage would be identical matrices), trigger Mode A
    # redundantly, and double memory usage for no benefit.

    local_anc_arr <- gf_windowed_pca(
      dosage_mat  = dosage_mat,   # 1KGP3 — PCA fitted here directly
      snp_chroms  = snp_chroms,
      k_pops      = K_POPS,
      window_size = 500L,
      step_size   = 250L,         # 50% overlap for smoother estimates
      scale_snps  = FALSE,
      n_cores     = 1L
    )
    cat("Windowed PCA: Mode B (study-only, no separate reference).\n")

  } else {

    # ── External study cohort — Mode A ────────────────────────────────────────
    # dosage_mat is an external cohort. Load 1KGP3 as reference so per-window
    # PCA axes are anchored to 1KGP3 population structure, not the study cohort.
    # Study samples are projected into each window's reference PC subspace.

    ref_raw    <- data.table::fread(qc$dosage_raw, nThread = 4L)
    ref_dosage <- as.matrix(ref_raw[, -(1:6)])
    ref_dosage[is.na(ref_dosage)] <- 0
    storage.mode(ref_dosage) <- "numeric"
    rm(ref_raw)

    local_anc_arr <- gf_windowed_pca(
      dosage_mat     = dosage_mat,   # external study cohort — projected in
      snp_chroms     = snp_chroms,
      ref_dosage_mat = ref_dosage,   # 1KGP3 — PCA fitted here per window
      k_pops         = K_POPS,
      window_size    = 500L,
      step_size      = 250L,
      scale_snps     = FALSE,
      n_cores        = 1L
    )
    cat("Windowed PCA: Mode A (reference-based, study cohort projected in).\n")
  }

  # Save per-window prcomp objects — required for scoring new individuals
  gf_save_window_pcas(local_anc_arr, file.path(OUTDIR, "window_pcas.rds"))
  cat("window_pcas.rds saved.\n")

} else if (LOCAL_ANC_METHOD == "flare") {

  # ── FLARE — external study cohort only ────────────────────────────────────
  # Validated above to require IS_1KG_COHORT = FALSE.
  # Convert QC'd pgen to VCF for use as FLARE reference.

  system2(PLINK2, c(
    "--pfile",        qc$pgen_prefix,
    "--recode",       "vcf", "bgz",
    "--real-ref-alleles",
    "--out",          file.path(OUTDIR, "plink", "kg3_ref"),
    "--silent"
  ))
  system2("tabix", c("-p", "vcf",
                     file.path(OUTDIR, "plink", "kg3_ref.vcf.gz")))
  KG3_REF_VCF <- file.path(OUTDIR, "plink", "kg3_ref.vcf.gz")

  flare_out <- gf_local_ancestry_flare(
    vcf_file   = STUDY_VCF,       # external study cohort (unphased OK)
    ref_vcf    = KG3_REF_VCF,     # 1KGP3 — independent reference
    sample_map = REF_PANEL,
    gmap_dir   = GMAP_DIR,
    outdir     = file.path(OUTDIR, "flare"),
    chroms     = CHROMS,
    k_pops     = K_POPS,
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
    global_anc        = anc$anc_probs,
    entropy_threshold = 0.3
  )
  # FLARE's K axes are named biological populations anchored to REF_PANEL.
  # No window_pcas.rds equivalent is needed. At scoring time, provide the same
  # REF_PANEL and KG3_REF_VCF to gf_local_ancestry_flare() and the coordinate
  # system is automatically consistent.

} else if (LOCAL_ANC_METHOD == "none") {
  cat("Global-only conditioning — local_anc will be NULL in gf_train().\n")
}

if (!is.null(local_anc_arr)) {
  saveRDS(local_anc_arr, file.path(OUTDIR, "local_anc_array.rds"))
  cat(sprintf("Local ancestry array: [%s]\n",
              paste(dim(local_anc_arr), collapse = " x ")))

  # Sanity check: row sums should equal 1 at every (sample, SNP)
  rs <- apply(local_anc_arr, c(1,2), sum)
  cat(sprintf("Row-sum check: min=%.4f  max=%.4f  (expect 1.0)\n",
              min(rs), max(rs)))

  # Per-chromosome summary
  la_summary <- gf_local_anc_summary(local_anc_arr, snp_chroms, anc$labels)
  write.csv(la_summary, file.path(OUTDIR, "local_anc_summary.csv"),
            row.names = FALSE)
}


# =============================================================================
# PHASE 7  Build and train the transformer
# =============================================================================
cat("[7/8] Building GenoFormer...\n")

use_local <- !is.null(local_anc_arr)
k_model   <- if (use_local) dim(local_anc_arr)[3] else K_POPS

model <- gf_build_model(
  n_snps        = M,                          # sequence length
  n_ld_blocks   = max(ld_blocks) + 1L,        # LD embedding table size
  d_model       = D_MODEL,
  n_heads       = N_HEADS,
  n_layers      = N_LAYERS,
  anc_dim       = ncol(anc$anc_probs),        # 5 superpopulations
  k_pops        = k_model,                    # K from local ancestry
  d_local       = D_LOCAL,
  use_local_anc = use_local,
  dropout       = 0.1,
  device        = "auto"                      # cuda > mps > cpu
)
print(model)

cat(sprintf("[7/8] Training — %d epochs, batch=%d, lr=%g...\n",
            EPOCHS, BATCH_SIZE, LR))

# Handle missing phenotype rows
valid <- !is.na(phenotype)
if (sum(valid) < N)
  cat(sprintf("Training on %d/%d samples with complete phenotype.\n",
              sum(valid), N))

# local_anc argument to gf_train:
#
#   When use_local = TRUE (windowed_pca or flare):
#     local_anc_arr[valid, , , drop=FALSE] — shape (N_valid x M x K)
#     At each training step, a batch of B rows is drawn: shape (B x M x K).
#     The LocalAncestryProjector MLP maps (B x M x K) → (B x M x d_local),
#     which is concatenated with the dosage and beta projections in each SNP
#     token embedding before attention. The model learns ancestry-stratified
#     LD and effect-size patterns simultaneously.
#
#     For IS_1KG_COHORT=TRUE with windowed_pca, K = K_POPS (pseudo-ancestry
#     dimensions from per-window PCA). For FLARE, K = number of populations
#     in REF_PANEL (named biological populations, e.g. AFR/AMR/EAS/EUR/SAS).
#     Both produce valid (B x M x K) tensors; only the interpretation of the
#     K dimension differs.
#
#   When use_local = FALSE (method="none", or local_anc_arr is NULL):
#     local_anc = NULL. gf_build_model was called with use_local_anc=FALSE so
#     no LocalAncestryProjector MLP exists. The model conditions only on global
#     ancestry via FiLM layers. k_pops and d_local parameters are ignored.

trained <- gf_train(
  model        = model,
  dosage_mat   = dosage_mat[valid, , drop=FALSE],
  anc_probs    = anc$anc_probs[valid, , drop=FALSE],
  pgs_weights  = pgs_weights,
  ld_blocks    = ld_blocks,
  chroms       = snp_chroms,
  phenotype    = phenotype[valid],
  pop_labels   = pop_int[valid],
  local_anc    = if (use_local) local_anc_arr[valid, , , drop=FALSE] else NULL,
  epochs       = EPOCHS,
  batch_size   = BATCH_SIZE,
  lr           = LR,
  weight_decay = WEIGHT_DECAY,
  anc_weight   = ANC_WEIGHT,
  calib_weight = CALIB_WEIGHT,
  val_split    = VAL_SPLIT,
  save_path    = file.path(OUTDIR, "genoformer_best.pt"),
  verbose      = TRUE
)
model <- trained$model

# Save model + architecture config
gf_save_model(
  model,
  path        = file.path(OUTDIR, "genoformer_best.pt"),
  save_config = TRUE    # writes genoformer_best_config.json
)
cat(sprintf("Model saved: genoformer_best.pt + genoformer_best_config.json\n"))

# Save training history
write.csv(trained$history,
          file.path(OUTDIR, "training_history.csv"),
          row.names = FALSE)

# Training loss plot
history <- trained$history
p_loss  <- ggplot(history, aes(x = epoch)) +
  geom_line(aes(y = train_loss, colour = "Train")) +
  geom_line(aes(y = val_loss,   colour = "Validation")) +
  scale_colour_manual(values = c(Train="#2E75B6", Validation="#E36C09"),
                      name = NULL) +
  labs(title = sprintf("GenoFormer training (%s)", LOCAL_ANC_METHOD),
       x = "Epoch", y = "Loss") +
  theme_minimal(base_size = 12)
ggplot2::ggsave(file.path(OUTDIR, "training_loss.png"), p_loss,
                width = 8, height = 5, dpi = 150)

cat(sprintf("Best validation loss: %.6f\n",
            min(history$val_loss, na.rm = TRUE)))


# =============================================================================
# PHASE 8  Compute training PRS and evaluate portability
# =============================================================================
cat("[8/8] Computing training PRS and evaluating portability...\n")

prs_train <- gf_predict(
  model          = model,
  dosage_mat     = dosage_mat,
  anc_probs      = anc$anc_probs,
  pgs_weights    = pgs_weights,
  ld_blocks      = ld_blocks,
  chroms         = snp_chroms,
  effect_alleles = effect_alleles,
  local_anc      = local_anc_arr,
  batch_size     = 128L
)
prs_train$sample_id  <- sample_ids
prs_train$population <- anc$labels

# Save training PRS scores
# Population means and SDs in this file are used to compute percentiles for
# new individuals at scoring time
saveRDS(prs_train, file.path(OUTDIR, "train_prs_scores.rds"))
write.csv(prs_train, file.path(OUTDIR, "train_prs_scores.csv"),
          row.names = FALSE)

# Correlation with phenotype per population
cat("\nCorrelation with phenotype by population:\n")
for (pop in c("AFR","AMR","EAS","EUR","SAS")) {
  mask <- prs_train$population == pop & !is.na(phenotype)
  if (sum(mask) < 3L) next
  r_t <- cor(prs_train$transformer_prs[mask], phenotype[mask], use="complete.obs")
  r_c <- cor(prs_train$classical_prs[mask],   phenotype[mask], use="complete.obs")
  cat(sprintf("  %-4s  transformer r=%.3f  classical r=%.3f\n", pop, r_t, r_c))
}

# Portability evaluation (requires phenotype for all samples)
pheno_eval <- phenotype
pheno_eval[is.na(pheno_eval)] <- mean(phenotype, na.rm=TRUE)  # impute for eval only

eval_res <- gf_evaluate(
  prs_tbl     = prs_train,
  phenotype   = pheno_eval,
  pop_labels  = anc$labels,
  ref_pop     = "EUR",
  n_bootstrap = 500L
)
cat("\nPortability metrics:\n")
print(eval_res$metrics)

ggplot2::ggsave(file.path(OUTDIR, "portability.png"),
                eval_res$portability_plot,  width = 8, height = 5, dpi = 150)
ggplot2::ggsave(file.path(OUTDIR, "distributions.png"),
                eval_res$distribution_plot, width = 8, height = 6, dpi = 150)


# =============================================================================
# SUMMARY
# =============================================================================

cat("\n════════════════════════════════════════════════\n")
cat(" Model generation complete\n")
cat("════════════════════════════════════════════════\n")
cat(sprintf(" Training cohort:    %s\n",
            if (IS_1KG_COHORT) "1KGP3 base model" else "External study cohort"))
cat(sprintf(" Output directory:   %s\n", OUTDIR))
cat(sprintf(" Samples trained:    %d\n", sum(valid)))
cat(sprintf(" SNPs (M):           %d\n", M))
cat(sprintf(" LD blocks:          %d\n", max(ld_blocks) + 1L))
cat(sprintf(" Model parameters:   %s\n", format(model$n_params, big.mark=",")))
cat(sprintf(" Local ancestry:     %s\n", LOCAL_ANC_METHOD))
cat(sprintf(" Device used:        %s\n", model$device))
cat(sprintf(" Best val loss:      %.6f\n", min(history$val_loss, na.rm=TRUE)))
cat("\n Files saved:\n")

artefacts <- list(
  "genoformer_best.pt"          = "Model weights",
  "genoformer_best_config.json" = "Architecture config",
  "ancestry.rds"                = "PCA + RF (global ancestry)",
  "window_pcas.rds"             = "Per-window PCAs (windowed_pca only)",
  "pgs_snps.txt"                = "PGS SNP list for new VCF extraction",
  "train_prs_scores.rds"        = "Training PRS (percentile normalisation)",
  "training_history.csv"        = "Loss per epoch",
  "training_loss.png"           = "Loss curve",
  "portability.png"             = "Cross-pop R² portability",
  "distributions.png"           = "PRS distributions by population"
)

for (fname in names(artefacts)) {
  full_path <- file.path(OUTDIR, fname)
  exists    <- file.exists(full_path)
  marker    <- if (exists) "✓" else "-"
  cat(sprintf("   %s %-36s %s\n", marker, fname, artefacts[[fname]]))
}

cat("\n genoformerR version:", as.character(packageVersion("genoformerR")), "\n")


# =============================================================================
# OPTIONAL: use gf_run_pipeline() for a shorter equivalent
# =============================================================================
# gf_run_pipeline() wraps all eight phases into one call. Note that v0.3.1
# manual steps above (ref_dosage_raw=, gf_save_window_pcas, gf_save_model)
# are not yet reflected in gf_run_pipeline; use the step-by-step script above
# for full v0.3.1 compatibility.
#
# result <- gf_run_pipeline(
#   dosage_raw       = qc$dosage_raw,
#   pgs_file         = PGS_FILE,
#   panel_file       = PANEL_FILE,
#   phenotype        = PHENOTYPE,
#   outdir           = OUTDIR,
#   local_anc_method = LOCAL_ANC_METHOD,   # "flare" | "windowed_pca" | "none"
#   vcf_file         = STUDY_VCF,          # FLARE only
#   ref_vcf          = KG3_REF_VCF,        # FLARE only
#   sample_map       = REF_PANEL,          # FLARE only
#   gmap_dir         = GMAP_DIR,           # FLARE only
#   flare_jar        = FLARE_JAR,          # FLARE only
#   k_pops           = K_POPS,
#   model_config     = list(d_model=D_MODEL, n_heads=N_HEADS, n_layers=N_LAYERS),
#   train_config     = list(epochs=EPOCHS, batch_size=BATCH_SIZE, lr=LR),
#   chroms           = CHROMS,
#   seed             = 42L
# )
