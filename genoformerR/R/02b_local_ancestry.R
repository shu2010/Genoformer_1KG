# genoformerR/R/02b_local_ancestry.R  (v0.3.0)
# Local ancestry inference — unified dispatcher supporting:
#   "flare"        — FLARE (Java, no phasing required, recommended)
#   "rfmix2"       — RFMix2 (legacy, requires SHAPEIT4 pre-phasing)
#   "windowed_pca" — pure R windowed PCA proxy (no external tools)
# All methods return the same (N × M × K) array consumed by gf_train / gf_predict.


# ══════════════════════════════════════════════════════════════════════════════
#  UNIFIED DISPATCHER
# ══════════════════════════════════════════════════════════════════════════════

#' Infer local ancestry using FLARE, RFMix2, or windowed PCA
#'
#' Unified dispatcher for per-SNP local ancestry inference. Recommended method
#' is \code{"flare"} — it accepts unphased genotypes directly, requires no
#' pre-phasing step, and is faster than RFMix2 at comparable accuracy.
#' \code{"windowed_pca"} requires no external binaries or reference panel and
#' is suitable for exploratory analyses or cohorts without a matched reference.
#'
#' @param method Character. One of \code{"flare"} (default), \code{"rfmix2"},
#'   or \code{"windowed_pca"}.
#' @param dosage_mat   Numeric matrix (N × M). Required for \code{"windowed_pca"};
#'   not used by FLARE/RFMix2 (those read VCFs directly).
#' @param vcf_file     Character. Path to study VCF (unphased OK for FLARE;
#'   must be phased for RFMix2). Required for \code{"flare"} and \code{"rfmix2"}.
#' @param ref_vcf      Character. Reference panel VCF. Required for FLARE and
#'   RFMix2.
#' @param sample_map   Character. Path to 2-column TSV (sample_id, population).
#'   Reference samples only. Required for FLARE and RFMix2.
#' @param gmap_dir     Character. Genetic map directory. Required for FLARE and
#'   RFMix2.
#' @param outdir       Character. Output directory. Default \code{"data/local_anc"}.
#' @param chroms       Integer vector. Chromosomes to process. Default 1:22.
#' @param k_pops       Integer. Expected number of reference populations K.
#'   Inferred from \code{sample_map} if NULL.
#' @param snp_ids      Character vector (M). SNP IDs in dosage matrix order.
#'   Required to align output to dosage SNP order.
#' @param snp_positions Integer vector (M). Physical positions. Required for
#'   alignment.
#' @param snp_chroms   Integer vector (M). Chromosome per SNP. Required for
#'   alignment.
#' @param sample_ids   Character vector (N). Sample IDs in dosage matrix order.
#' @param global_anc   Numeric matrix (N × K). Global ancestry from
#'   \code{\link{gf_ancestry}}, used as fallback for non-admixed samples.
#' @param ...          Additional arguments passed to the method-specific
#'   function (\code{\link{gf_local_ancestry_flare}},
#'   \code{\link{gf_local_ancestry_rfmix2}}, or
#'   \code{\link{gf_windowed_pca}}).
#'
#' @return A numeric array of shape \code{(N, M, K)} with dimnames set on all
#'   axes. K = number of reference populations. Rows sum to 1 at each SNP.
#' @export
#'
#' @examples
#' \dontrun{
#' # FLARE (recommended — no phasing needed)
#' local_anc <- gf_local_ancestry(
#'   method       = "flare",
#'   vcf_file     = "study.vcf.gz",
#'   ref_vcf      = "kg3_ref.vcf.gz",
#'   sample_map   = "sample_map.txt",
#'   gmap_dir     = "genetic_maps",
#'   snp_ids      = colnames(dosage_mat),
#'   snp_positions = snp_pos,
#'   snp_chroms   = snp_chr,
#'   sample_ids   = rownames(dosage_mat),
#'   global_anc   = anc$anc_probs
#' )
#'
#' # Windowed PCA (no external tools required)
#' local_anc <- gf_local_ancestry(
#'   method      = "windowed_pca",
#'   dosage_mat  = dosage_mat,
#'   snp_chroms  = snp_chr,
#'   k_pops      = 5L
#' )
#' }
gf_local_ancestry <- function(method       = c("flare", "rfmix2", "windowed_pca"),
                               dosage_mat   = NULL,
                               vcf_file     = NULL,
                               ref_vcf      = NULL,
                               sample_map   = NULL,
                               gmap_dir     = NULL,
                               outdir       = "data/local_anc",
                               chroms       = 1:22,
                               k_pops       = NULL,
                               snp_ids      = NULL,
                               snp_positions = NULL,
                               snp_chroms   = NULL,
                               sample_ids   = NULL,
                               global_anc   = NULL,
                               ...) {
  method <- match.arg(method)
  cli::cli_h1("Local Ancestry Inference — {toupper(method)}")

  switch(method,
    flare = {
      .require_args(list(vcf_file = vcf_file, ref_vcf = ref_vcf,
                         sample_map = sample_map, gmap_dir = gmap_dir,
                         snp_ids = snp_ids, snp_positions = snp_positions,
                         snp_chroms = snp_chroms, sample_ids = sample_ids),
                    method = "flare")
      raw <- gf_local_ancestry_flare(
        vcf_file = vcf_file, ref_vcf = ref_vcf,
        sample_map = sample_map, gmap_dir = gmap_dir,
        outdir = outdir, chroms = chroms, k_pops = k_pops, ...
      )
      gf_align_local_anc(raw, snp_ids, snp_positions, snp_chroms,
                          sample_ids, global_anc)
    },
    rfmix2 = {
      .require_args(list(vcf_file = vcf_file, ref_vcf = ref_vcf,
                         sample_map = sample_map, gmap_dir = gmap_dir,
                         snp_ids = snp_ids, snp_positions = snp_positions,
                         snp_chroms = snp_chroms, sample_ids = sample_ids),
                    method = "rfmix2")
      raw <- gf_local_ancestry_rfmix2(
        vcf_file = vcf_file, ref_vcf = ref_vcf,
        sample_map = sample_map, gmap_dir = gmap_dir,
        outdir = outdir, chroms = chroms, k_pops = k_pops, ...
      )
      gf_align_local_anc(raw, snp_ids, snp_positions, snp_chroms,
                          sample_ids, global_anc)
    },
    windowed_pca = {
      .require_args(list(dosage_mat = dosage_mat, snp_chroms = snp_chroms),
                    method = "windowed_pca")
      k_pops_wp <- if (is.null(k_pops)) 5L else as.integer(k_pops)
      gf_windowed_pca(dosage_mat = dosage_mat, snp_chroms = snp_chroms,
                      k_pops = k_pops_wp, ...)
    }
  )
}


# ══════════════════════════════════════════════════════════════════════════════
#  METHOD 1: FLARE  (recommended)
# ══════════════════════════════════════════════════════════════════════════════

#' Run FLARE local ancestry inference
#'
#' FLARE (Fast Local Ancestry REsolution) infers per-SNP ancestry posteriors
#' directly from unphased genotype data — no pre-phasing with SHAPEIT4
#' required. It is substantially faster than RFMix2 and is the recommended
#' method for biobank-scale datasets.
#'
#' @details
#' FLARE requires Java ≥ 11. Download the jar from
#' \url{https://github.com/browning-lab/flare}.
#'
#' Output files written by FLARE:
#' \itemize{
#'   \item \code{<prefix>.anc.vcf.gz} — per-haplotype ancestry calls (FORMAT:
#'     \code{ANS} field: two space-separated ancestry integers per sample)
#'   \item \code{<prefix>.log} — run log
#' }
#'
#' @param vcf_file   Character. Path to study VCF (unphased accepted).
#' @param ref_vcf    Character. Path to reference panel VCF.
#' @param sample_map Character. 2-column headerless TSV: sample_id, population.
#' @param gmap_dir   Character. Directory of genetic maps (PLINK-format
#'   \code{.gmap} or \code{chr*.txt}).
#' @param outdir     Character. Output directory. Default \code{"data/flare"}.
#' @param chroms     Integer vector. Chromosomes to run. Default 1:22.
#' @param k_pops     Integer or NULL. Expected K; validated against sample_map.
#' @param flare_jar  Character. Path to \code{flare.jar}. Default \code{"flare.jar"}.
#' @param java       Character. Path to java binary. Default \code{"java"}.
#' @param min_maf    Numeric. Minimum MAF filter passed to FLARE. Default 0.005.
#' @param nthreads   Integer. Parallel threads. Default 4.
#' @param em_its     Integer. EM iterations. Default 3.
#'
#' @return Named list of per-chromosome FLARE output prefixes.
#' @export
gf_local_ancestry_flare <- function(vcf_file,
                                     ref_vcf,
                                     sample_map,
                                     gmap_dir,
                                     outdir    = "data/flare",
                                     chroms    = 1:22,
                                     k_pops    = NULL,
                                     flare_jar = "flare.jar",
                                     java      = "java",
                                     min_maf   = 0.005,
                                     nthreads  = 4L,
                                     em_its    = 3L) {

  .check_java(java)
  .check_jar(flare_jar, "https://github.com/browning-lab/flare")
  fs::dir_create(outdir)

  # Validate reference population count from sample_map
  smap    <- data.table::fread(sample_map, header = FALSE, nThread = 2L)
  k_found <- length(unique(smap[[2]]))
  if (!is.null(k_pops) && k_found != k_pops) {
    cli::cli_warn("sample_map has {k_found} populations but k_pops={k_pops}. Using {k_found}.")
  }
  k_pops <- k_found
  cli::cli_alert_info("FLARE: K={k_pops} reference populations")
  cli::cli_alert_info("Note: FLARE accepts unphased VCF — no pre-phasing required")

  out_prefixes <- list()
  pb <- cli::cli_progress_bar("FLARE chromosomes", total = length(chroms))

  for (chr in chroms) {
    chr_c   <- as.character(chr)
    out_pfx <- file.path(outdir, glue::glue("chr{chr}"))

    # Locate genetic map for this chromosome
    gmap <- .find_gmap(gmap_dir, chr)
    if (is.na(gmap)) {
      cli::cli_warn("No genetic map for chr{chr}. Skipping.")
      cli::cli_progress_update(); next
    }

    # FLARE command-line arguments
    # gt=     study genotype VCF (unphased accepted)
    # ref=    reference panel VCF
    # map=    genetic map file
    # ref-panel= sample map file
    # out=    output prefix
    args <- c(
      "-jar", flare_jar,
      glue::glue("gt={vcf_file}"),
      glue::glue("ref={ref_vcf}"),
      glue::glue("map={gmap}"),
      glue::glue("ref-panel={sample_map}"),
      glue::glue("out={out_pfx}"),
      glue::glue("chrom={chr}"),
      glue::glue("nthreads={nthreads}"),
      glue::glue("min-maf={min_maf}"),
      glue::glue("em-its={em_its}")
    )

    ret <- system2(java, args, stdout = FALSE, stderr = FALSE)
    if (ret != 0L) {
      cli::cli_warn("FLARE exited with code {ret} on chr{chr}. Skipping.")
    } else {
      out_prefixes[[chr_c]] <- out_pfx
      cli::cli_alert_success(
        "chr{chr} \u2192 {.path {out_pfx}.anc.vcf.gz}"
      )
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  cli::cli_alert_success("FLARE complete: {length(out_prefixes)} chromosome(s).")
  invisible(out_prefixes)
}


#' Parse FLARE output into an aligned (N × M × K) local ancestry array
#'
#' Reads FLARE \code{.anc.vcf.gz} files and aligns the per-SNP, per-haplotype
#' ancestry calls to the dosage SNP order. The \code{ANS} FORMAT field
#' contains two integers per sample (one per haplotype); these are averaged to
#' produce diploid posteriors. One-hot encoding is applied (FLARE gives hard
#' calls by default; soft posteriors are produced if \code{em-its > 0}).
#'
#' @param flare_prefixes Named list of FLARE output prefixes per chromosome
#'   (output of \code{\link{gf_local_ancestry_flare}}).
#' @param snp_ids        Character vector (M). SNP IDs in dosage matrix order.
#' @param snp_positions  Integer vector (M). Physical positions.
#' @param snp_chroms     Integer vector (M). Chromosome per SNP.
#' @param sample_ids     Character vector (N). Sample IDs in dosage matrix order.
#' @param global_anc     Numeric matrix (N × K). Global ancestry fallback for
#'   non-admixed samples.
#' @param entropy_threshold Numeric. Admixture entropy threshold below which
#'   global ancestry is used instead of local. Default 0.3.
#' @param n_cores        Integer. Parallel cores for parsing. Default 2.
#'
#' @return Numeric array (N × M × K).
#' @export
gf_parse_flare <- function(flare_prefixes,
                            snp_ids,
                            snp_positions,
                            snp_chroms,
                            sample_ids,
                            global_anc          = NULL,
                            entropy_threshold   = 0.3,
                            n_cores             = 2L) {

  cli::cli_h1("Parsing FLARE Output")

  # Infer population names from first available .anc.vcf.gz
  first_pfx   <- flare_prefixes[[1]]
  anc_vcf     <- paste0(first_pfx, ".anc.vcf.gz")
  if (!file.exists(anc_vcf))
    cli::cli_abort("FLARE output not found: {.path {anc_vcf}}")

  hdr <- .read_vcf_header(anc_vcf)
  # FLARE embeds population labels in the VCF header:
  # ##ancestry=AFR,AMR,EAS,EUR,SAS
  anc_line  <- hdr[grepl("^##ancestry=", hdr)]
  if (length(anc_line) == 0) {
    cli::cli_warn("Cannot find ##ancestry header. Using generic Pop1..PopK labels.")
    pop_names <- NULL
  } else {
    pop_names <- strsplit(sub("^##ancestry=", "", anc_line[1]), ",")[[1]]
  }
  K <- length(pop_names)
  cli::cli_alert_info("K = {K} populations: {paste(pop_names, collapse=', ')}")

  N <- length(sample_ids)
  M <- length(snp_ids)

  # Pre-fill with global ancestry as baseline
  out <- .init_array(N, M, K, pop_names, sample_ids, snp_ids, global_anc)

  admixed_mask <- .admixed_samples(global_anc, entropy_threshold, N)
  cli::cli_alert_info(
    "{sum(admixed_mask)} admixed / {sum(!admixed_mask)} non-admixed samples"
  )

  pb <- cli::cli_progress_bar("Parsing FLARE chromosomes",
                               total = length(flare_prefixes))
  for (chr_c in names(flare_prefixes)) {
    chr     <- as.integer(chr_c)
    anc_vcf <- paste0(flare_prefixes[[chr_c]], ".anc.vcf.gz")
    if (!file.exists(anc_vcf)) {
      cli::cli_warn("Missing: {.path {anc_vcf}}"); cli::cli_progress_update(); next
    }

    chr_mask  <- snp_chroms == chr
    if (!any(chr_mask)) { cli::cli_progress_update(); next }

    vcf <- .read_flare_anc_vcf(anc_vcf, sample_ids, pop_names, K)
    # vcf: list with $pos (integer), $anc_array (n_vcf_snps × N × K) diploid

    # Nearest-SNP alignment
    dosage_pos  <- snp_positions[chr_mask]
    dosage_idx  <- which(chr_mask)
    nearest     <- .nearest_idx(dosage_pos, vcf$pos)

    for (n_i in which(admixed_mask)) {
      out[n_i, dosage_idx, ] <- vcf$anc_array[nearest, n_i, ]
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  out <- .normalise_rows(out)
  cli::cli_alert_success(
    "FLARE array ready: [{N} \u00d7 {M} \u00d7 {K}]"
  )
  out
}


# ══════════════════════════════════════════════════════════════════════════════
#  METHOD 2: RFMix2  (legacy — requires pre-phasing)
# ══════════════════════════════════════════════════════════════════════════════

#' Phase haplotypes with SHAPEIT4
#'
#' Required pre-processing step when using \code{method = "rfmix2"}. Not needed
#' for FLARE or windowed PCA.
#'
#' @param pgen_prefix Character. PLINK2 pgen prefix.
#' @param outdir      Character. Output directory for phased VCFs.
#' @param gmap_dir    Character. Genetic map directory.
#' @param chroms      Integer vector. Default 1:22.
#' @param shapeit4    Character. Path to shapeit4 binary. Default \code{"shapeit4"}.
#' @param scaffold    Character or NULL. Scaffold VCF for reference-based phasing.
#' @param ncores      Integer. Threads per chromosome. Default 4.
#' @param pbwt_depth  Integer. PBWT depth. Default 4.
#'
#' @return Named character vector of phased VCF paths (one per chromosome).
#' @export
gf_phase <- function(pgen_prefix,
                     outdir     = "data/phased",
                     gmap_dir,
                     chroms     = 1:22,
                     shapeit4   = "shapeit4",
                     scaffold   = NULL,
                     ncores     = 4L,
                     pbwt_depth = 4L) {

  .check_binary(shapeit4, "https://odelaneau.github.io/shapeit4/")
  fs::dir_create(outdir)
  cli::cli_h1("Haplotype Phasing — SHAPEIT4")
  phased_vcfs <- character(0)

  for (chr in chroms) {
    gmap <- .find_gmap(gmap_dir, chr)
    if (is.na(gmap)) {
      cli::cli_warn("No genetic map for chr{chr}. Skipping."); next
    }
    out_vcf <- file.path(outdir, glue::glue("chr{chr}_phased.vcf.gz"))
    tmp_vcf <- tempfile(fileext = ".vcf.gz")
    .run_plink2("plink2", "--pfile", pgen_prefix, "--chr", chr,
                "--recode", "vcf", "bgz",
                "--out", sub("\\.vcf\\.gz$", "", tmp_vcf), "--silent")

    args <- c("--input", tmp_vcf, "--map", gmap, "--output", out_vcf,
              "--region", glue::glue("chr{chr}"),
              "--thread", ncores, "--pbwt-depth", pbwt_depth)
    if (!is.null(scaffold))
      args <- c(args, "--scaffold", scaffold,
                "--scaffold-region", glue::glue("chr{chr}"))

    ret <- system2(shapeit4, args, stdout = FALSE, stderr = FALSE)
    if (ret != 0L) {
      cli::cli_warn("SHAPEIT4 exited with code {ret} on chr{chr}.")
    } else {
      system2("tabix", c("-p", "vcf", out_vcf), stdout = FALSE, stderr = FALSE)
      phased_vcfs[as.character(chr)] <- out_vcf
      cli::cli_alert_success("chr{chr} phased")
    }
    unlink(tmp_vcf)
  }
  cli::cli_alert_success("Phasing complete: {length(phased_vcfs)} chromosome(s).")
  invisible(phased_vcfs)
}


#' Run RFMix2 local ancestry inference (legacy)
#'
#' Requires phased VCFs from \code{\link{gf_phase}}. Prefer
#' \code{\link{gf_local_ancestry_flare}} for new analyses.
#'
#' @param vcf_file     Character. Path to phased study VCF (or a named vector
#'   of per-chromosome phased VCFs; the \code{chrom} argument selects one).
#' @param ref_vcf      Character. Reference panel VCF.
#' @param sample_map   Character. 2-column TSV (sample_id, population).
#' @param gmap_dir     Character. Genetic map directory.
#' @param outdir       Character. Output directory. Default \code{"data/rfmix2"}.
#' @param chroms       Integer vector. Chromosomes to process. Default 1:22.
#' @param k_pops       Integer or NULL. Expected K.
#' @param rfmix2       Character. Path to rfmix2 binary. Default \code{"rfmix2"}.
#' @param em_iterations Integer. EM iterations. Default 1.
#'
#' @return Named list of RFMix2 output prefixes per chromosome.
#' @export
gf_local_ancestry_rfmix2 <- function(vcf_file,
                                      ref_vcf,
                                      sample_map,
                                      gmap_dir,
                                      outdir        = "data/rfmix2",
                                      chroms        = 1:22,
                                      k_pops        = NULL,
                                      rfmix2        = "rfmix2",
                                      em_iterations = 1L) {

  .check_binary(rfmix2, "https://github.com/slowkoni/rfmix")
  fs::dir_create(outdir)
  cli::cli_alert_info("RFMix2 (legacy): requires pre-phased VCFs")

  if (is.null(k_pops)) {
    smap   <- data.table::fread(sample_map, header = FALSE, nThread = 2L)
    k_pops <- length(unique(smap[[2]]))
  }
  cli::cli_alert_info("K={k_pops} reference populations")

  out_prefixes <- list()
  pb <- cli::cli_progress_bar("RFMix2 chromosomes", total = length(chroms))

  for (chr in chroms) {
    chr_c  <- as.character(chr)
    gmap   <- .find_gmap(gmap_dir, chr)
    if (is.na(gmap)) {
      cli::cli_warn("No genetic map for chr{chr}. Skipping.")
      cli::cli_progress_update(); next
    }

    # Support per-chromosome VCF vector OR single multi-chromosome VCF
    study_vcf <- if (length(vcf_file) > 1 && chr_c %in% names(vcf_file)) {
      vcf_file[chr_c]
    } else {
      vcf_file[[1]]
    }

    out_pfx <- file.path(outdir, glue::glue("chr{chr}"))
    args <- c(
      "-f", study_vcf, "-r", ref_vcf, "-m", sample_map, "-g", gmap,
      "-o", out_pfx,
      "--chromosome",    glue::glue("chr{chr}"),
      "--n-iterations",  em_iterations,
      "--reanalyze-reference"
    )
    ret <- system2(rfmix2, args, stdout = FALSE, stderr = FALSE)
    if (ret != 0L) {
      cli::cli_warn("RFMix2 exited with code {ret} on chr{chr}.")
    } else {
      out_prefixes[[chr_c]] <- out_pfx
      cli::cli_alert_success("chr{chr} \u2192 {.path {out_pfx}.fb.tsv}")
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  invisible(out_prefixes)
}


# ══════════════════════════════════════════════════════════════════════════════
#  METHOD 3: WINDOWED PCA  (pure R, no external tools)
# ══════════════════════════════════════════════════════════════════════════════

#' Windowed PCA proxy for local ancestry
#'
#' Computes principal components within sliding genomic windows and converts
#' them to pseudo-ancestry proportions via softmax normalisation.
#'
#' @details
#' ## Reference-based mode (\code{ref_dosage_mat} or \code{ref_window_pcas} provided)
#'
#' PCA is fitted on the reference samples (e.g. 1KGP3) per window. Study
#' samples are then projected into that fixed reference PC space using the saved
#' rotation matrix. This ensures the k pseudo-ancestry axes have the same
#' meaning at training time and at scoring time, regardless of the size or
#' composition of the study cohort. Without this, the PC axes at scoring time
#' are arbitrary rotations with no relationship to the axes seen during training.
#'
#' Per-window \code{prcomp} objects produced during training should be saved
#' with \code{\link{gf_save_window_pcas}} and loaded at scoring time with
#' \code{\link{gf_load_window_pcas}}, then passed as \code{ref_window_pcas=}.
#'
#' ## Study-only mode (\code{ref_dosage_mat = NULL}, \code{ref_window_pcas = NULL})
#'
#' PCA is fitted on \code{dosage_mat} directly. Only valid when \code{dosage_mat}
#' already is the 1KGP3 training cohort (N large, stable per-window PCs). The
#' resulting window PCA objects must still be saved for consistent scoring later.
#' Not suitable for single-sample or small-cohort scoring.
#'
#' @param dosage_mat      Numeric matrix (N x M). Study dosage values.
#' @param snp_chroms      Integer vector (M). Chromosome per SNP.
#' @param ref_dosage_mat  Numeric matrix (N_ref x M) or NULL. Reference dosage
#'   (e.g. 1KGP3). When provided, PCA is fitted on these samples and
#'   \code{dosage_mat} is projected in. Recommended for any external cohort.
#' @param ref_window_pcas Named list of saved \code{prcomp} objects from a
#'   prior training run (from \code{\link{gf_load_window_pcas}}). Takes
#'   precedence over \code{ref_dosage_mat}. Use this when scoring new
#'   individuals after training.
#' @param k_pops          Integer. Number of pseudo-ancestry dimensions. Default 5.
#' @param window_size     Integer. SNPs per window. Default 500.
#' @param step_size       Integer. Step between window starts. Defaults to
#'   \code{window_size} (non-overlapping).
#' @param min_var         Numeric. Minimum variance explained by PC1. Default 0.01.
#' @param scale_snps      Logical. Scale SNPs to unit variance. Default FALSE.
#' @param n_cores         Integer. Parallel cores. Default 1.
#'
#' @return Numeric array (N x M x k_pops) with an attribute \code{"window_pcas"}
#'   containing per-window \code{prcomp} objects. Save with
#'   \code{\link{gf_save_window_pcas}}.
#' @export
#'
#' @examples
#' \dontrun{
#' # Training on 1KGP3 — study-only mode (dosage_mat IS 1KGP3)
#' la_train <- gf_windowed_pca(kg3_dosage, snp_chroms)
#' gf_save_window_pcas(la_train, "models/window_pcas.rds")
#'
#' # Scoring a new cohort — project into saved training subspace
#' ref_pcas  <- gf_load_window_pcas("models/window_pcas.rds")
#' la_new    <- gf_windowed_pca(new_dosage, snp_chroms,
#'                               ref_window_pcas = ref_pcas)
#'
#' # Training on a non-1KGP3 cohort — reference-based mode
#' la_train  <- gf_windowed_pca(study_dosage, snp_chroms,
#'                               ref_dosage_mat = kg3_dosage)
#' gf_save_window_pcas(la_train, "models/window_pcas.rds")
#' }
gf_windowed_pca <- function(dosage_mat,
                             snp_chroms,
                             ref_dosage_mat  = NULL,
                             ref_window_pcas = NULL,
                             k_pops          = 5L,
                             window_size     = 500L,
                             step_size       = NULL,
                             min_var         = 0.01,
                             scale_snps      = FALSE,
                             n_cores         = 1L) {

  N <- nrow(dosage_mat)
  M <- ncol(dosage_mat)
  step <- if (is.null(step_size)) window_size else as.integer(step_size)

  use_saved_pcas <- !is.null(ref_window_pcas)
  use_ref_mat    <- !is.null(ref_dosage_mat) && !use_saved_pcas

  if (use_saved_pcas) {
    # k_pops inferred from saved PCAs
    k_pops <- ncol(ref_window_pcas[[1]]$rotation)
    cli::cli_alert_info(
      "Windowed PCA — projection mode: {length(ref_window_pcas)} saved window PCAs, K={k_pops}."
    )
  } else if (use_ref_mat) {
    N_ref  <- nrow(ref_dosage_mat)
    k_pops <- min(as.integer(k_pops), N_ref - 1L)
    cli::cli_alert_info(
      "Windowed PCA — reference mode: fitting on {N_ref} reference samples, projecting {N} study sample(s), K={k_pops}."
    )
  } else {
    k_pops <- min(as.integer(k_pops), N - 1L)
    if (N < 30L)
      cli::cli_warn(c(
        "N={N} is very small for study-only windowed PCA.",
        "i" = "PC axes will be unstable and inconsistent with training.",
        "i" = "Provide ref_dosage_mat= (e.g. 1KGP3) or ref_window_pcas= to project into a stable reference subspace."
      ))
    cli::cli_alert_info(
      "Windowed PCA — study-only mode: fitting on {N} samples, K={k_pops}."
    )
  }

  pop_names   <- paste0("wPC", seq_len(k_pops))
  out         <- array(0, dim = c(N, M, k_pops),
                       dimnames = list(rownames(dosage_mat),
                                       colnames(dosage_mat),
                                       pop_names))
  new_win_pcas <- list()

  chrs <- sort(unique(snp_chroms))
  cli::cli_alert_info(
    "window={window_size} SNPs, step={step}, K={k_pops}"
  )
  pb <- cli::cli_progress_bar("Windowed PCA chromosomes", total = length(chrs))

  for (chr in chrs) {
    chr_mask <- snp_chroms == chr
    chr_idx  <- which(chr_mask)
    n_snps   <- length(chr_idx)

    if (n_snps < k_pops + 1L) {
      cli::cli_warn("chr{chr}: only {n_snps} SNPs < k_pops+1. Skipping.")
      cli::cli_progress_update(); next
    }

    starts <- seq(1L, n_snps - window_size + 1L, by = step)
    if (length(starts) == 0L) starts <- 1L

    run_window <- function(s) {
      w_end   <- min(s + window_size - 1L, n_snps)
      w_snps  <- chr_idx[s:w_end]
      win_key <- paste0(chr, ":", s)

      if (use_saved_pcas) {
        # ── Mode C: project into pre-saved per-window PCAs ───────────────────
        pca <- ref_window_pcas[[win_key]]
        if (is.null(pca)) return(NULL)
        # align study columns to the SNPs the saved PCA was fitted on
        ref_snp_names <- rownames(pca$rotation)
        study_snps    <- colnames(dosage_mat)[w_snps]
        shared        <- intersect(study_snps, ref_snp_names)
        if (length(shared) < k_pops) return(NULL)
        scores <- stats::predict(pca, newdata = dosage_mat[, shared, drop = FALSE])[, seq_len(k_pops), drop = FALSE]
        pca_out <- NULL   # don't re-store saved PCAs

      } else if (use_ref_mat) {
        # ── Mode A: fit on reference, project study ───────────────────────────
        ref_w  <- ref_dosage_mat[, w_snps, drop = FALSE]
        keep   <- apply(ref_w, 2, stats::var) > 0
        if (sum(keep) < k_pops) return(NULL)
        ref_w  <- ref_w[, keep, drop = FALSE]
        pca    <- tryCatch(
          stats::prcomp(ref_w, center = TRUE, scale. = scale_snps, rank. = k_pops),
          error = function(e) NULL
        )
        if (is.null(pca)) return(NULL)
        if ((pca$sdev[1]^2 / sum(pca$sdev^2)) < min_var) return(NULL)
        scores  <- stats::predict(pca, newdata = dosage_mat[, w_snps[keep], drop = FALSE])[, seq_len(k_pops), drop = FALSE]
        pca_out <- pca

      } else {
        # ── Mode B: study-only PCA ────────────────────────────────────────────
        geno_w <- dosage_mat[, w_snps, drop = FALSE]
        keep   <- apply(geno_w, 2, stats::var) > 0
        if (sum(keep) < k_pops) return(NULL)
        geno_w <- geno_w[, keep, drop = FALSE]
        pca    <- tryCatch(
          stats::prcomp(geno_w, center = TRUE, scale. = scale_snps, rank. = k_pops),
          error = function(e) NULL
        )
        if (is.null(pca)) return(NULL)
        if ((pca$sdev[1]^2 / sum(pca$sdev^2)) < min_var) return(NULL)
        scores  <- pca$x[, seq_len(k_pops), drop = FALSE]
        pca_out <- pca
      }

      # Softmax normalisation: shift to non-negative, then exp/rowSums
      scores <- scores - apply(scores, 1, min)
      exp_s  <- exp(scores)
      soft_s <- exp_s / rowSums(exp_s)
      list(snps = w_snps, soft = soft_s, key = win_key, pca = pca_out)
    }

    if (n_cores > 1L) {
      cl      <- parallel::makeCluster(min(n_cores, length(starts)))
      results <- parallel::parLapply(cl, starts, run_window)
      parallel::stopCluster(cl)
    } else {
      results <- lapply(starts, run_window)
    }

    counts <- numeric(M)
    for (res in results) {
      if (is.null(res)) next
      for (k in seq_len(k_pops)) {
        out[, res$snps, k] <- out[, res$snps, k] + res$soft[, k]
      }
      counts[res$snps] <- counts[res$snps] + 1L
      if (!is.null(res$pca)) new_win_pcas[[res$key]] <- res$pca
    }

    covered <- counts > 0
    for (m in which(chr_mask)) {
      if (covered[m]) {
        out[, m, ] <- out[, m, ] / counts[m]
      } else {
        cov_idx <- chr_idx[covered[chr_idx]]
        if (length(cov_idx) == 0L) {
          out[, m, ] <- 1 / k_pops
        } else {
          near <- cov_idx[which.min(abs(cov_idx - m))]
          out[, m, ] <- out[, near, ]
        }
      }
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  out <- .normalise_rows(out)

  if (length(new_win_pcas) > 0) {
    attr(out, "window_pcas") <- new_win_pcas
    if (!use_saved_pcas)
      cli::cli_alert_info(
        "Save {length(new_win_pcas)} window PCAs with gf_save_window_pcas() for consistent scoring of new individuals."
      )
  }

  cli::cli_alert_success(
    "Windowed PCA array ready: [{N} \u00d7 {M} \u00d7 {k_pops}]"
  )
  out
}


#' Save per-window PCA objects for future projection
#'
#' Extracts the \code{"window_pcas"} attribute set by
#' \code{\link{gf_windowed_pca}} and writes it to an RDS file. These objects
#' must be saved to ensure that new individuals can be projected into the same
#' local ancestry subspace that was used during training.
#'
#' @param local_anc_arr Array returned by \code{\link{gf_windowed_pca}}.
#' @param path          Character. Output \code{.rds} path.
#' @export
gf_save_window_pcas <- function(local_anc_arr,
                                 path = "models/window_pcas.rds") {
  pcas <- attr(local_anc_arr, "window_pcas")
  if (is.null(pcas) || length(pcas) == 0)
    cli::cli_abort(c(
      "No window_pcas attribute on this array.",
      "i" = "Run gf_windowed_pca() without ref_window_pcas= to generate and attach window PCA objects."
    ))
  saveRDS(pcas, path)
  cli::cli_alert_success("Saved {length(pcas)} window PCAs to {.path {path}}")
  invisible(path)
}


#' Load saved per-window PCA objects
#'
#' Loads the RDS file written by \code{\link{gf_save_window_pcas}}.
#' Pass the result as \code{ref_window_pcas=} to \code{\link{gf_windowed_pca}}
#' when computing local ancestry for new individuals.
#'
#' @param path Character. Path to the \code{.rds} file.
#' @return Named list of \code{prcomp} objects keyed by \code{"chr:window_start"}.
#' @export
gf_load_window_pcas <- function(path = "models/window_pcas.rds") {
  if (!file.exists(path))
    cli::cli_abort("Window PCA file not found: {.path {path}}")
  pcas <- readRDS(path)
  cli::cli_alert_success("Loaded {length(pcas)} window PCAs from {.path {path}}")
  pcas
}


# ══════════════════════════════════════════════════════════════════════════════
#  SHARED ALIGNMENT  (FLARE + RFMix2)
# ══════════════════════════════════════════════════════════════════════════════

#' Align RFMix2 .fb.tsv posteriors to dosage SNP order
#'
#' Used by \code{gf_local_ancestry(method = "rfmix2")} after running
#' \code{\link{gf_local_ancestry_rfmix2}}. Handles nearest-SNP interpolation
#' and entropy-based fallback to global ancestry for non-admixed samples.
#'
#' @param rfmix_prefixes Named list of RFMix2 output prefixes.
#' @param snp_ids        Character vector (M).
#' @param snp_positions  Integer vector (M).
#' @param snp_chroms     Integer vector (M).
#' @param sample_ids     Character vector (N).
#' @param global_anc     Numeric matrix (N × K). Fallback ancestry.
#' @param entropy_threshold Numeric. Default 0.3.
#' @param n_cores        Integer. Default 2.
#'
#' @return Numeric array (N × M × K).
#' @export
gf_align_rfmix2 <- function(rfmix_prefixes,
                              snp_ids,
                              snp_positions,
                              snp_chroms,
                              sample_ids,
                              global_anc,
                              entropy_threshold = 0.3,
                              n_cores           = 2L) {

  cli::cli_h1("Aligning RFMix2 Posteriors")
  first_pfx <- rfmix_prefixes[[1]]
  fb_file   <- paste0(first_pfx, ".fb.tsv")
  if (!file.exists(fb_file))
    cli::cli_abort("RFMix2 output not found: {.path {fb_file}}")

  # Infer K and population names from header
  fb_hdr    <- data.table::fread(fb_file, nrows = 0L, nThread = 2L)
  hdr       <- names(fb_hdr)
  pop_names <- unique(sub(".*:::hap[01]:::", "", hdr[grepl(":::", hdr)]))
  K         <- length(pop_names)

  N  <- length(sample_ids)
  M  <- length(snp_ids)
  out <- .init_array(N, M, K, pop_names, sample_ids, snp_ids, global_anc)
  admixed_mask <- .admixed_samples(global_anc, entropy_threshold, N)

  pb <- cli::cli_progress_bar("Aligning RFMix2 chromosomes",
                               total = length(rfmix_prefixes))
  for (chr_c in names(rfmix_prefixes)) {
    chr      <- as.integer(chr_c)
    fb_file  <- paste0(rfmix_prefixes[[chr_c]], ".fb.tsv")
    if (!file.exists(fb_file)) {
      cli::cli_warn("Missing: {.path {fb_file}}"); cli::cli_progress_update(); next
    }
    chr_mask   <- snp_chroms == chr
    if (!any(chr_mask)) { cli::cli_progress_update(); next }

    fb         <- data.table::fread(fb_file, nThread = n_cores)
    pos_col    <- if ("#pos" %in% names(fb)) "#pos" else "pos"
    rfmix_pos  <- fb[[pos_col]]
    nearest    <- .nearest_idx(snp_positions[chr_mask], rfmix_pos)
    dosage_idx <- which(chr_mask)

    for (n_i in which(admixed_mask)) {
      sid <- sample_ids[n_i]
      for (k_i in seq_len(K)) {
        pop      <- pop_names[k_i]
        col_hap0 <- paste0(sid, ":::hap0:::", pop)
        col_hap1 <- paste0(sid, ":::hap1:::", pop)
        if (!col_hap0 %in% names(fb)) next
        dip <- (as.numeric(fb[[col_hap0]]) + as.numeric(fb[[col_hap1]])) / 2
        out[n_i, dosage_idx, k_i] <- dip[nearest]
      }
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()
  .normalise_rows(out)
}


#' Generic alignment helper — dispatches to FLARE or RFMix2 parser
#' @noRd
gf_align_local_anc <- function(raw_prefixes, snp_ids, snp_positions,
                                snp_chroms, sample_ids, global_anc,
                                entropy_threshold = 0.3, n_cores = 2L) {
  # Detect format by checking which output file exists
  first <- raw_prefixes[[1]]
  if (file.exists(paste0(first, ".anc.vcf.gz"))) {
    gf_parse_flare(raw_prefixes, snp_ids, snp_positions, snp_chroms,
                   sample_ids, global_anc, entropy_threshold, n_cores)
  } else if (file.exists(paste0(first, ".fb.tsv"))) {
    gf_align_rfmix2(raw_prefixes, snp_ids, snp_positions, snp_chroms,
                    sample_ids, global_anc, entropy_threshold, n_cores)
  } else {
    cli::cli_abort("Cannot detect output format at {.path {first}}")
  }
}


#' Summary statistics for a local ancestry array
#'
#' Reports mean ancestry per population per chromosome, and flags samples with
#' suspiciously uniform local ancestry (possible inference failure).
#'
#' @param local_anc_arr Array (N × M × K).
#' @param snp_chroms    Integer vector (M).
#' @param pop_labels    Character vector (N). Known superpopulation.
#'
#' @return Data frame with columns: labelled_pop, chr, and one column per K.
#' @export
gf_local_anc_summary <- function(local_anc_arr, snp_chroms, pop_labels) {
  K     <- dim(local_anc_arr)[3]
  pops  <- dimnames(local_anc_arr)[[3]]
  chrs  <- sort(unique(snp_chroms))
  rows  <- list()
  for (pop in unique(pop_labels)) {
    pop_mask <- pop_labels == pop
    for (chr in chrs) {
      chr_mask <- snp_chroms == chr
      if (!any(chr_mask)) next
      mean_probs <- colMeans(
        apply(local_anc_arr[pop_mask, chr_mask, , drop = FALSE], c(1, 3), mean)
      )
      row        <- as.list(round(mean_probs, 4))
      names(row) <- pops
      rows[[length(rows) + 1]] <- c(list(labelled_pop = pop, chr = chr), row)
    }
  }
  do.call(rbind, lapply(rows, as.data.frame))
}


# ══════════════════════════════════════════════════════════════════════════════
#  INTERNAL HELPERS
# ══════════════════════════════════════════════════════════════════════════════

.check_binary <- function(bin, url) {
  if (Sys.which(bin) == "" && !file.exists(bin))
    cli::cli_abort(c("Binary {.val {bin}} not found.", "i" = "Install: {url}"))
}

.check_java <- function(java = "java") {
  if (Sys.which(java) == "" && !file.exists(java))
    cli::cli_abort("Java not found. Install Java >= 11: https://adoptium.net")
  # Version check
  out <- tryCatch(
    system2(java, "-version", stdout = TRUE, stderr = TRUE),
    error = function(e) ""
  )
  maj <- suppressWarnings(as.integer(
    sub(".*version \"([0-9]+).*\".*", "\\1", paste(out, collapse = " "))
  ))
  if (!is.na(maj) && maj < 11L)
    cli::cli_abort("FLARE requires Java >= 11. Found Java {maj}.")
}

.check_jar <- function(jar, url) {
  if (!file.exists(jar))
    cli::cli_abort(c("JAR not found: {.path {jar}}", "i" = "Download: {url}"))
}

.find_gmap <- function(gmap_dir, chr) {
  hits <- fs::dir_ls(gmap_dir,
                     regexp = glue::glue("chr{chr}[._].*\\.(txt|gmap)(\\.gz)?$"))
  if (length(hits) == 0L) NA_character_ else hits[[1]]
}

.nearest_idx <- function(query_pos, ref_pos) {
  vapply(query_pos, function(qp) which.min(abs(ref_pos - qp)), integer(1L))
}

.init_array <- function(N, M, K, pop_names, sample_ids, snp_ids, global_anc) {
  out <- array(NA_real_, dim = c(N, M, K),
               dimnames = list(sample_ids, snp_ids, pop_names))
  if (!is.null(global_anc)) {
    k_match <- match(pop_names, colnames(global_anc))
    valid_k <- !is.na(k_match)
    for (m in seq_len(M))
      out[, m, valid_k] <- global_anc[, k_match[valid_k], drop = FALSE]
  }
  out
}

.admixed_samples <- function(global_anc, threshold, N) {
  if (is.null(global_anc)) return(rep(TRUE, N))
  entropy <- -.rowSums(
    global_anc * log(global_anc + 1e-9),
    nrow(global_anc), ncol(global_anc)
  )
  entropy >= threshold
}

.normalise_rows <- function(arr) {
  row_sums <- apply(arr, c(1, 2), sum)
  for (k in seq_len(dim(arr)[3]))
    arr[, , k] <- arr[, , k] / pmax(row_sums, 1e-8)
  arr
}

.require_args <- function(args_list, method) {
  missing <- names(Filter(is.null, args_list))
  if (length(missing) > 0)
    cli::cli_abort(
      "method='{method}' requires: {paste(missing, collapse=', ')}"
    )
}

# Read VCF header lines (starting with #)
.read_vcf_header <- function(vcf_gz) {
  con   <- gzcon(file(vcf_gz, "rb"))
  lines <- character(0)
  repeat {
    ln <- readLines(con, n = 1L, warn = FALSE)
    if (length(ln) == 0L || !startsWith(ln, "#")) break
    lines <- c(lines, ln)
  }
  close(con)
  lines
}

# Parse FLARE .anc.vcf.gz: returns list(pos, anc_array)
# anc_array shape: (n_vcf_snps, N, K) diploid averaged
.read_flare_anc_vcf <- function(anc_vcf, sample_ids, pop_names, K) {
  # Read full VCF (skip ## comment lines, keep # CHROM header)
  raw      <- data.table::fread(
    cmd     = paste("zcat", anc_vcf),
    skip    = "#CHROM",
    nThread = 2L
  )
  pos      <- raw[["POS"]]
  fmt_idx  <- which(names(raw) == "FORMAT")
  samp_cols <- names(raw)[(fmt_idx + 1):ncol(raw)]

  # Match VCF sample columns to requested sample_ids
  N        <- length(sample_ids)
  n_vcf    <- nrow(raw)
  anc_arr  <- array(0, dim = c(n_vcf, N, K))

  # FLARE FORMAT field: ANS (ancestry per haplotype, space-separated integers)
  # E.g. "0 1" means hap0=pop0, hap1=pop1
  ans_field_idx <- function(format_str) {
    flds <- strsplit(format_str, ":")[[1]]
    which(flds == "ANS")
  }
  fmt_ex  <- raw[["FORMAT"]][1]
  ans_idx <- ans_field_idx(fmt_ex)

  for (n_i in seq_along(sample_ids)) {
    sid <- sample_ids[n_i]
    if (!sid %in% samp_cols) next
    gts <- raw[[sid]]

    # Extract ANS field from GT string "GT:ANS:..." -> "0 1" etc.
    ans_vals <- vapply(strsplit(gts, ":"), function(x) {
      if (length(x) >= ans_idx) x[[ans_idx]] else NA_character_
    }, character(1L))

    # Parse two haplotype integers
    hap_mat <- do.call(rbind, lapply(strsplit(ans_vals, " "), function(x) {
      suppressWarnings(as.integer(x[1:2]))
    }))

    # Convert to diploid proportions across K populations
    for (h in 1:2) {
      hap <- hap_mat[, h]  # integer ancestry index (0-based)
      for (s in seq_len(n_vcf)) {
        if (!is.na(hap[s]) && hap[s] < K)
          anc_arr[s, n_i, hap[s] + 1L] <- anc_arr[s, n_i, hap[s] + 1L] + 0.5
      }
    }
  }
  list(pos = pos, anc_array = anc_arr)
}
