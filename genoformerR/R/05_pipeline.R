# genoformerR/R/05_pipeline.R  (v0.3.0)
# End-to-end wrapper — supports local_anc_method: "flare", "rfmix2",
# "windowed_pca", or "none"

#' Run the complete GenoFormer ancestry-aware PRS pipeline
#'
#' Orchestrates all steps from raw PLINK2 dosage files to evaluated PRS scores.
#' v0.3 replaces the separate \code{run_phasing} / \code{run_local_anc} flags
#' with a single \code{local_anc_method} argument.
#'
#' Supported local ancestry methods:
#' \describe{
#'   \item{\code{"flare"}}{Recommended. Java FLARE — no pre-phasing required,
#'     fast, biobank-scale. Requires \code{vcf_file}, \code{ref_vcf},
#'     \code{sample_map}, \code{gmap_dir}.}
#'   \item{\code{"rfmix2"}}{Legacy. Requires SHAPEIT4 pre-phasing and
#'     \code{vcf_file} (phased), \code{ref_vcf}, \code{sample_map},
#'     \code{gmap_dir}.}
#'   \item{\code{"windowed_pca"}}{Pure R — no external tools or reference panel.
#'     Slower for very large M, suitable for exploratory analyses.}
#'   \item{\code{"none"}}{Global-only conditioning (v0.1 behaviour).}
#' }
#'
#' @param dosage_raw        Character. Path to PLINK2 \code{.raw} dosage file.
#' @param pgs_file          Character. Path to PGS Catalog scoring file.
#' @param panel_file        Character. Path to 1KGP3 sample panel file.
#' @param phenotype         Numeric vector or path to phenotype file (2-col: ID, pheno).
#' @param outdir            Character. Output directory. Default \code{"genoformer_out"}.
#' @param local_anc_method  Character. One of \code{"flare"} (default),
#'   \code{"rfmix2"}, \code{"windowed_pca"}, or \code{"none"}.
#' @param vcf_file          Character or NULL. Study VCF (unphased for FLARE;
#'   phased for RFMix2). Required for \code{"flare"} and \code{"rfmix2"}.
#' @param ref_vcf           Character or NULL. Reference panel VCF. Required for
#'   \code{"flare"} and \code{"rfmix2"}.
#' @param sample_map        Character or NULL. 2-col TSV: sample_id, population.
#'   Required for \code{"flare"} and \code{"rfmix2"}.
#' @param gmap_dir          Character or NULL. Genetic map directory. Required for
#'   \code{"flare"} and \code{"rfmix2"}.
#' @param k_pops            Integer. Reference population count K. Default 5.
#'   For \code{"windowed_pca"}, sets the number of pseudo-ancestry dimensions.
#' @param flare_jar         Character. Path to \code{flare.jar}. Default
#'   \code{"flare.jar"}. Only used when \code{local_anc_method = "flare"}.
#' @param run_phasing       Logical. Run SHAPEIT4 before RFMix2. Only relevant
#'   when \code{local_anc_method = "rfmix2"}. Default FALSE.
#' @param pgen_prefix       Character or NULL. PLINK2 pgen prefix for phasing.
#'   Required when \code{run_phasing = TRUE}.
#' @param model_config      Named list. Passed to \code{\link{gf_build_model}}.
#' @param train_config      Named list. Passed to \code{\link{gf_train}}.
#' @param ancestry_method   Character. Global ancestry method: \code{"r"} or
#'   \code{"python"}. Default \code{"r"}.
#' @param device            Character. \code{"auto"}, \code{"cuda"}, or
#'   \code{"cpu"}. Default \code{"auto"}.
#' @param chroms            Integer vector. Chromosomes to process. Default 1:22.
#' @param seed              Integer. Random seed. Default 42.
#'
#' @return Named list: \code{ancestry}, \code{local_anc_arr}, \code{model},
#'   \code{history}, \code{prs}, \code{evaluation}, \code{paths}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Global-only (no external tools)
#' gf_run_pipeline(dosage_raw="...", pgs_file="...", panel_file="...",
#'                 phenotype=pheno, local_anc_method="none")
#'
#' # FLARE (recommended for admixed cohorts)
#' gf_run_pipeline(dosage_raw="...", pgs_file="...", panel_file="...",
#'                 phenotype=pheno, local_anc_method="flare",
#'                 vcf_file="study.vcf.gz", ref_vcf="ref.vcf.gz",
#'                 sample_map="ref_panel.txt", gmap_dir="gmaps/",
#'                 flare_jar="flare.jar")
#'
#' # Windowed PCA (no external dependencies)
#' gf_run_pipeline(dosage_raw="...", pgs_file="...", panel_file="...",
#'                 phenotype=pheno, local_anc_method="windowed_pca",
#'                 k_pops=5L)
#' }
gf_run_pipeline <- function(dosage_raw,
                             pgs_file,
                             panel_file,
                             phenotype,
                             outdir           = "genoformer_out",
                             local_anc_method = c("flare", "windowed_pca",
                                                   "rfmix2", "none"),
                             vcf_file         = NULL,
                             ref_vcf          = NULL,
                             sample_map       = NULL,
                             gmap_dir         = NULL,
                             k_pops           = 5L,
                             flare_jar        = "flare.jar",
                             run_phasing      = FALSE,
                             pgen_prefix      = NULL,
                             model_config     = list(),
                             train_config     = list(),
                             ancestry_method  = "r",
                             device           = "auto",
                             chroms           = 1:22,
                             seed             = 42L) {

  local_anc_method <- match.arg(local_anc_method)
  set.seed(seed)
  fs::dir_create(outdir)
  t0 <- proc.time()

  # ── Mode banner ─────────────────────────────────────────────────────────────
  mode_label <- c(
    flare        = "FLARE (no phasing required)",
    rfmix2       = "RFMix2 (legacy, requires phasing)",
    windowed_pca = "Windowed PCA (pure R, no reference panel)",
    none         = "Global-only (no local ancestry)"
  )
  cli::cli_h1("GenoFormer v0.3 Pipeline")
  cli::cli_alert_info("Local ancestry: {.strong {mode_label[local_anc_method]}}")

  # ── Validate method-specific inputs ─────────────────────────────────────────
  if (local_anc_method %in% c("flare", "rfmix2")) {
    missing_args <- names(Filter(is.null,
      list(vcf_file = vcf_file, ref_vcf = ref_vcf,
           sample_map = sample_map, gmap_dir = gmap_dir)))
    if (length(missing_args) > 0)
      cli::cli_abort(
        "local_anc_method='{local_anc_method}' requires: {paste(missing_args, collapse=', ')}"
      )
  }

  # ── 0. Python backend ────────────────────────────────────────────────────────
  gf_init()

  # ── 1. Load dosage + PGS ─────────────────────────────────────────────────────
  cli::cli_alert_info("[1/6] Loading dosage matrix...")
  raw        <- data.table::fread(dosage_raw, nThread = 4L)
  sample_ids <- raw$IID
  dosage_mat <- as.matrix(raw[, -(1:6)])
  dosage_mat[is.na(dosage_mat)] <- 0L

  pgs       <- data.table::fread(pgs_file, skip = "#", nThread = 2L)
  names(pgs) <- tolower(names(pgs))
  w_col    <- intersect(c("effect_weight", "beta", "or"),     names(pgs))[1]
  ea_col   <- intersect(c("effect_allele", "a1"),             names(pgs))[1]
  id_col   <- intersect(c("rsid", "variant_id", "id", "snp"), names(pgs))[1]

  snp_ids        <- colnames(dosage_mat)
  pgs_matched    <- pgs[match(snp_ids, pgs[[id_col]]), ]
  pgs_weights    <- as.numeric(pgs_matched[[w_col]])
  effect_alleles <- as.character(pgs_matched[[ea_col]])
  chroms_snp     <- .snp_chroms(snp_ids)
  snp_positions  <- .snp_positions(snp_ids)
  ld_blocks      <- .assign_ld_blocks(snp_ids, chroms_snp)

  if (is.character(phenotype) && length(phenotype) == 1L && file.exists(phenotype)) {
    ph_df     <- data.table::fread(phenotype, nThread = 2L)
    phenotype <- as.numeric(ph_df[[2]][match(sample_ids, ph_df[[1]])])
  }
  stopifnot(length(phenotype) == length(sample_ids))

  # ── 2. Global ancestry ────────────────────────────────────────────────────────
  cli::cli_alert_info("[2/6] Inferring global ancestry...")
  anc <- gf_ancestry(dosage_mat = dosage_mat, panel_file = panel_file,
                     method = ancestry_method)
  saveRDS(anc, file.path(outdir, "ancestry.rds"))

  # ── 3. Local ancestry ─────────────────────────────────────────────────────────
  cli::cli_alert_info("[3/6] Local ancestry — {local_anc_method}...")
  local_anc_arr <- NULL

  if (local_anc_method == "flare") {
    raw_prefixes <- gf_local_ancestry_flare(
      vcf_file   = vcf_file,
      ref_vcf    = ref_vcf,
      sample_map = sample_map,
      gmap_dir   = gmap_dir,
      outdir     = file.path(outdir, "flare"),
      chroms     = chroms,
      k_pops     = k_pops,
      flare_jar  = flare_jar
    )
    local_anc_arr <- gf_parse_flare(
      flare_prefixes = raw_prefixes,
      snp_ids        = snp_ids,
      snp_positions  = snp_positions,
      snp_chroms     = chroms_snp,
      sample_ids     = sample_ids,
      global_anc     = anc$anc_probs
    )

  } else if (local_anc_method == "rfmix2") {
    # Optional phasing step
    if (run_phasing) {
      if (is.null(pgen_prefix))
        cli::cli_abort("pgen_prefix required when run_phasing=TRUE")
      phased_vcfs <- gf_phase(
        pgen_prefix = pgen_prefix,
        outdir      = file.path(outdir, "phased"),
        gmap_dir    = gmap_dir,
        chroms      = chroms
      )
      study_vcf <- phased_vcfs
    } else {
      study_vcf <- vcf_file
    }
    raw_prefixes <- gf_local_ancestry_rfmix2(
      vcf_file   = study_vcf,
      ref_vcf    = ref_vcf,
      sample_map = sample_map,
      gmap_dir   = gmap_dir,
      outdir     = file.path(outdir, "rfmix2"),
      chroms     = chroms,
      k_pops     = k_pops
    )
    local_anc_arr <- gf_align_rfmix2(
      rfmix_prefixes = raw_prefixes,
      snp_ids        = snp_ids,
      snp_positions  = snp_positions,
      snp_chroms     = chroms_snp,
      sample_ids     = sample_ids,
      global_anc     = anc$anc_probs
    )

  } else if (local_anc_method == "windowed_pca") {
    local_anc_arr <- gf_windowed_pca(
      dosage_mat = dosage_mat,
      snp_chroms = chroms_snp,
      k_pops     = k_pops
    )

  } else {
    cli::cli_alert_info("Local ancestry skipped — global-only conditioning.")
  }

  if (!is.null(local_anc_arr)) {
    saveRDS(local_anc_arr, file.path(outdir, "local_anc_array.rds"))
    dims <- paste(dim(local_anc_arr), collapse = " x ")
    cli::cli_alert_success("Local ancestry array: [{dims}]")
  }

  # ── 4. Build model ────────────────────────────────────────────────────────────
  cli::cli_alert_info("[4/6] Building GenoFormer...")
  use_local <- !is.null(local_anc_arr)
  # k_pops for model: infer from array if available
  k_model <- if (use_local) dim(local_anc_arr)[3] else k_pops
  model <- do.call(gf_build_model, c(
    list(n_snps        = ncol(dosage_mat),
         device        = device,
         use_local_anc = use_local,
         k_pops        = k_model),
    model_config
  ))

  # ── 5. Train ──────────────────────────────────────────────────────────────────
  cli::cli_alert_info("[5/6] Training GenoFormer...")
  pop_int <- as.integer(factor(anc$labels,
                               levels = c("AFR","AMR","EAS","EUR","SAS"))) - 1L
  trained <- do.call(gf_train, c(
    list(model        = model,
         dosage_mat   = dosage_mat,
         anc_probs    = anc$anc_probs,
         pgs_weights  = pgs_weights,
         ld_blocks    = ld_blocks,
         chroms       = chroms_snp,
         phenotype    = phenotype,
         pop_labels   = pop_int,
         local_anc    = local_anc_arr,
         save_path    = file.path(outdir, "genoformer_best.pt")),
    train_config
  ))
  model <- trained$model
  write.csv(trained$history, file.path(outdir, "training_history.csv"),
            row.names = FALSE)

  # ── 6. PRS + evaluation ───────────────────────────────────────────────────────
  cli::cli_alert_info("[6/6] Computing and evaluating PRS...")
  prs <- gf_predict(
    model          = model,
    dosage_mat     = dosage_mat,
    anc_probs      = anc$anc_probs,
    pgs_weights    = pgs_weights,
    ld_blocks      = ld_blocks,
    chroms         = chroms_snp,
    effect_alleles = effect_alleles,
    local_anc      = local_anc_arr
  )
  prs$sample_id  <- sample_ids
  prs$population <- anc$labels
  write.csv(prs, file.path(outdir, "prs_scores.csv"), row.names = FALSE)

  eval_res <- gf_evaluate(prs, phenotype, anc$labels)
  ggplot2::ggsave(file.path(outdir, "portability.png"),
                  eval_res$portability_plot,  width = 8, height = 5, dpi = 150)
  ggplot2::ggsave(file.path(outdir, "distributions.png"),
                  eval_res$distribution_plot, width = 8, height = 6, dpi = 150)

  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  cli::cli_alert_success(
    "Pipeline complete in {elapsed} min. Outputs: {.path {outdir}}"
  )

  list(
    ancestry      = anc,
    local_anc_arr = local_anc_arr,
    model         = model,
    history       = trained$history,
    prs           = prs,
    evaluation    = eval_res,
    paths         = list(
      outdir    = outdir,
      prs_csv   = file.path(outdir, "prs_scores.csv"),
      model_pt  = file.path(outdir, "genoformer_best.pt"),
      history   = file.path(outdir, "training_history.csv")
    )
  )
}


# ── Internal SNP metadata helpers ──────────────────────────────────────────────

.snp_chroms <- function(snp_ids) {
  parsed <- suppressWarnings(
    as.integer(sub("^chr([0-9XY]+)[_:].*", "\\1", snp_ids))
  )
  if (all(is.na(parsed))) return(rep_len(1:22, length(snp_ids)))
  ifelse(is.na(parsed), 1L, parsed)
}

.snp_positions <- function(snp_ids) {
  parsed <- suppressWarnings(
    as.integer(sub("^chr[0-9XY]+[_:]([0-9]+)[_:].*", "\\1", snp_ids))
  )
  if (all(is.na(parsed))) return(seq_along(snp_ids) * 1000L)
  ifelse(is.na(parsed), seq_along(snp_ids) * 1000L, parsed)
}

.assign_ld_blocks <- function(snp_ids, chroms, block_size = 200L) {
  block_within_chr <- stats::ave(
    seq_along(snp_ids), chroms,
    FUN = function(x) ceiling(seq_along(x) / block_size)
  )
  as.integer(factor(paste(chroms, block_within_chr, sep = "_"))) - 1L
}
