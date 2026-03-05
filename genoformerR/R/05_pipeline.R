# genoformerR/R/05_pipeline.R  (v0.2.0)
# End-to-end convenience wrapper with dual-level ancestry conditioning

#' Run the complete GenoFormer ancestry-aware PRS pipeline
#'
#' Orchestrates all steps from raw PLINK2 dosage files to evaluated PRS. In
#' v0.2, optionally runs SHAPEIT4 phasing and RFMix2 local ancestry inference
#' for dual-level (global + local) ancestry conditioning.
#'
#' @param dosage_raw      Character. Path to PLINK2 \code{.raw} dosage file.
#' @param pgs_file        Character. Path to PGS Catalog scoring file.
#' @param panel_file      Character. Path to 1KGP3 sample panel file.
#' @param phenotype       Numeric vector or path to phenotype file.
#' @param outdir          Character. Output directory.
#' @param run_phasing     Logical. Run SHAPEIT4 haplotype phasing. Default FALSE.
#' @param run_local_anc   Logical. Run RFMix2 local ancestry inference. Default FALSE.
#' @param gmap_dir        Character or NULL. Genetic map directory.
#' @param ref_vcf         Character or NULL. Reference panel VCF for RFMix2.
#' @param rfmix_sample_map Character or NULL. 2-col TSV: sample_id, population.
#' @param phased_vcf_dir  Character or NULL. Pre-phased VCFs directory.
#' @param k_pops          Integer. RFMix2 reference population count. Default 5.
#' @param model_config    Named list. Passed to \code{\link{gf_build_model}}.
#' @param train_config    Named list. Passed to \code{\link{gf_train}}.
#' @param ancestry_method Character. \code{"r"} or \code{"python"}.
#' @param device          Character. \code{"auto"}, \code{"cuda"}, or \code{"cpu"}.
#' @param chroms          Integer vector. Chromosomes to process. Default 1:22.
#' @param seed            Integer. Random seed. Default 42.
#'
#' @return Named list: ancestry, local_anc_arr, model, history, prs, evaluation, paths.
#' @export
gf_run_pipeline <- function(dosage_raw,
                             pgs_file,
                             panel_file,
                             phenotype,
                             outdir           = "genoformer_out",
                             run_phasing      = FALSE,
                             run_local_anc    = FALSE,
                             gmap_dir         = NULL,
                             ref_vcf          = NULL,
                             rfmix_sample_map = NULL,
                             phased_vcf_dir   = NULL,
                             k_pops           = 5L,
                             model_config     = list(),
                             train_config     = list(),
                             ancestry_method  = "r",
                             device           = "auto",
                             chroms           = 1:22,
                             seed             = 42L) {

  set.seed(seed)
  fs::dir_create(outdir)
  t0 <- proc.time()

  cli::cli_h1("GenoFormer v0.2 End-to-End Pipeline")
  mode_str <- if (run_local_anc) "Dual-level (global + local ancestry)" else "Global-only ancestry"
  cli::cli_alert_info("Mode: {.strong {mode_str}}")

  # 0. Init ──────────────────────────────────────────────────────────────────
  gf_init()

  # 1. Load dosage + PGS ─────────────────────────────────────────────────────
  cli::cli_alert_info("[1/7] Loading dosage matrix...")
  raw        <- data.table::fread(dosage_raw, nThread = 4L)
  sample_ids <- raw$IID
  dosage_mat <- as.matrix(raw[, -(1:6)])
  dosage_mat[is.na(dosage_mat)] <- 0L

  pgs      <- data.table::fread(pgs_file, skip = "#", nThread = 2L)
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

  if (is.character(phenotype) && length(phenotype) == 1 && file.exists(phenotype)) {
    ph_df     <- data.table::fread(phenotype, nThread = 2L)
    phenotype <- as.numeric(ph_df[[2]][match(sample_ids, ph_df[[1]])])
  }
  stopifnot(length(phenotype) == length(sample_ids))

  # 2. Global ancestry ───────────────────────────────────────────────────────
  cli::cli_alert_info("[2/7] Inferring global ancestry...")
  anc <- gf_ancestry(dosage_mat = dosage_mat, panel_file = panel_file,
                     method = ancestry_method)
  saveRDS(anc, file.path(outdir, "ancestry.rds"))

  # 3. Phasing ───────────────────────────────────────────────────────────────
  phased_vcfs <- NULL
  if (run_phasing && run_local_anc) {
    if (is.null(gmap_dir)) cli::cli_abort("gmap_dir required when run_phasing=TRUE")
    cli::cli_alert_info("[3/7] Phasing haplotypes (SHAPEIT4)...")
    phased_vcfs <- gf_phase(
      pgen_prefix = sub("\\.raw$", "", dosage_raw),
      outdir      = file.path(outdir, "phased"),
      gmap_dir    = gmap_dir,
      chroms      = chroms
    )
  } else if (!is.null(phased_vcf_dir) && run_local_anc) {
    phased_vcfs <- fs::dir_ls(phased_vcf_dir, regexp = "chr[0-9]+_phased\\.vcf\\.gz$")
    names(phased_vcfs) <- sub(".*chr([0-9]+)_phased.*", "\\1", phased_vcfs)
    cli::cli_alert_info("[3/7] Pre-phased VCFs: {length(phased_vcfs)} chromosome(s)")
  } else {
    cli::cli_alert_info("[3/7] Phasing: skipped")
  }

  # 4. Local ancestry ────────────────────────────────────────────────────────
  local_anc_arr <- NULL
  if (run_local_anc) {
    if (is.null(phased_vcfs))
      cli::cli_abort("Phased VCFs required. Set run_phasing=TRUE or phased_vcf_dir.")
    if (is.null(ref_vcf) || is.null(rfmix_sample_map))
      cli::cli_abort("ref_vcf and rfmix_sample_map required when run_local_anc=TRUE")

    cli::cli_alert_info("[4/7] Running RFMix2...")
    rfmix_out <- gf_local_ancestry(
      phased_vcfs  = phased_vcfs,
      ref_vcf      = ref_vcf,
      sample_map   = rfmix_sample_map,
      outdir       = file.path(outdir, "rfmix"),
      gmap_dir     = gmap_dir,
      n_ref_pops   = k_pops,
      chroms       = chroms
    )
    cli::cli_alert_info("[4b/7] Aligning local ancestry posteriors...")
    local_anc_arr <- gf_load_local_anc(
      rfmix_prefixes = rfmix_out,
      snp_ids        = snp_ids,
      snp_positions  = snp_positions,
      snp_chroms     = chroms_snp,
      sample_ids     = sample_ids,
      global_anc     = anc$anc_probs
    )
    saveRDS(local_anc_arr, file.path(outdir, "local_anc_array.rds"))
  } else {
    cli::cli_alert_info("[4/7] Local ancestry: skipped")
  }

  # 5. Build model ───────────────────────────────────────────────────────────
  cli::cli_alert_info("[5/7] Building GenoFormer...")
  model <- do.call(gf_build_model, c(
    list(n_snps = ncol(dosage_mat), device = device,
         use_local_anc = run_local_anc, k_pops = k_pops),
    model_config
  ))

  # 6. Train ─────────────────────────────────────────────────────────────────
  cli::cli_alert_info("[6/7] Training GenoFormer...")
  pop_int <- as.integer(factor(anc$labels,
                               levels = c("AFR","AMR","EAS","EUR","SAS"))) - 1L
  trained <- do.call(gf_train, c(
    list(model       = model,
         dosage_mat  = dosage_mat,
         anc_probs   = anc$anc_probs,
         pgs_weights = pgs_weights,
         ld_blocks   = ld_blocks,
         chroms      = chroms_snp,
         phenotype   = phenotype,
         pop_labels  = pop_int,
         local_anc   = local_anc_arr,
         save_path   = file.path(outdir, "genoformer_best.pt")),
    train_config
  ))
  model <- trained$model
  write.csv(trained$history, file.path(outdir, "training_history.csv"), row.names = FALSE)

  # 7. PRS + evaluation ──────────────────────────────────────────────────────
  cli::cli_alert_info("[7/7] Computing and evaluating PRS...")
  prs <- gf_predict(
    model = model, dosage_mat = dosage_mat, anc_probs = anc$anc_probs,
    pgs_weights = pgs_weights, ld_blocks = ld_blocks, chroms = chroms_snp,
    effect_alleles = effect_alleles, local_anc = local_anc_arr
  )
  prs$sample_id  <- sample_ids
  prs$population <- anc$labels
  write.csv(prs, file.path(outdir, "prs_scores.csv"), row.names = FALSE)

  eval_res <- gf_evaluate(prs, phenotype, anc$labels)
  ggplot2::ggsave(file.path(outdir, "portability.png"),
                  eval_res$portability_plot, width = 8, height = 5, dpi = 150)
  ggplot2::ggsave(file.path(outdir, "distributions.png"),
                  eval_res$distribution_plot, width = 8, height = 6, dpi = 150)

  elapsed <- round((proc.time() - t0)["elapsed"] / 60, 1)
  cli::cli_alert_success("Pipeline complete in {elapsed} min. Outputs: {.path {outdir}}")

  list(ancestry = anc, local_anc_arr = local_anc_arr,
       model = model, history = trained$history,
       prs = prs, evaluation = eval_res,
       paths = list(outdir   = outdir,
                    prs_csv  = file.path(outdir, "prs_scores.csv"),
                    model_pt = file.path(outdir, "genoformer_best.pt"),
                    history  = file.path(outdir, "training_history.csv")))
}


# ─── Internal helpers ──────────────────────────────────────────────────────────

.snp_chroms <- function(snp_ids) {
  parsed <- suppressWarnings(as.integer(
    sub("^chr([0-9XY]+)[_:].*", "\\1", snp_ids)
  ))
  if (all(is.na(parsed))) return(rep_len(1:22, length(snp_ids)))
  ifelse(is.na(parsed), 1L, parsed)
}

# Parse physical position from chr:pos:ref:alt or chr_pos_ref_alt SNP IDs
.snp_positions <- function(snp_ids) {
  parsed <- suppressWarnings(as.integer(
    sub("^chr[0-9XY]+[_:]([0-9]+)[_:].*", "\\1", snp_ids)
  ))
  if (all(is.na(parsed))) {
    # Fall back to sequential positions (for tests / demo data)
    return(seq_along(snp_ids) * 1000L)
  }
  ifelse(is.na(parsed), seq_along(snp_ids) * 1000L, parsed)
}

.assign_ld_blocks <- function(snp_ids, chroms, block_size = 200L) {
  block_within_chr <- ave(seq_along(snp_ids), chroms,
                          FUN = function(x) ceiling(seq_along(x) / block_size))
  as.integer(factor(paste(chroms, block_within_chr, sep = "_"))) - 1L
}
