# genoformerR/R/02b_local_ancestry.R
# Haplotype phasing (SHAPEIT4) and local ancestry inference (RFMix2)
# Output: per-SNP local ancestry posterior array (N_samples × N_snps × K_pops)

#' Phase haplotypes with SHAPEIT4
#'
#' Runs SHAPEIT4 per chromosome on a PLINK2 pgen dataset and writes phased
#' VCFs. Phased haplotypes are required input for RFMix2 local ancestry
#' inference. A genetic map and (optionally) a scaffold of known phase from
#' a reference panel are required for accurate phasing.
#'
#' @param pgen_prefix Character. Prefix of PLINK2 pgen files (no extension).
#' @param outdir      Character. Output directory for phased VCFs.
#' @param gmap_dir    Character. Directory containing per-chromosome genetic
#'   maps in SHAPEIT4 format (e.g., from the Bherer et al. 2017 maps).
#' @param chroms      Integer vector. Chromosomes to phase. Default 1:22.
#' @param shapeit4    Character. Path to shapeit4 binary. Default \code{"shapeit4"}.
#' @param scaffold    Character or NULL. Path to scaffold VCF for reference-
#'   based phasing (e.g., 1KGP3 phased haplotypes). Greatly improves accuracy
#'   for admixed samples.
#' @param ncores      Integer. Threads per chromosome. Default 4.
#' @param pbwt_depth  Integer. PBWT depth (higher = more accurate, slower).
#'   Default 4.
#'
#' @return Invisibly returns a named character vector of phased VCF paths,
#'   one per chromosome.
#' @export
#'
#' @examples
#' \dontrun{
#' phased <- gf_phase(
#'   pgen_prefix = "data/plink/kg3_qc",
#'   outdir      = "data/phased",
#'   gmap_dir    = "data/genetic_maps",
#'   chroms      = 21:22
#' )
#' }
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
    cli::cli_alert_info("Phasing chr{chr}...")
    gmap <- fs::dir_ls(gmap_dir,
                       regexp = glue::glue("chr{chr}[._].*\\.txt(\\.gz)?$"))[1]
    if (is.na(gmap)) {
      cli::cli_warn("No genetic map found for chr{chr} in {gmap_dir}. Skipping.")
      next
    }
    out_vcf <- file.path(outdir, glue::glue("chr{chr}_phased.vcf.gz"))

    # Convert pgen → VCF for this chromosome
    tmp_vcf <- tempfile(fileext = ".vcf.gz")
    .run_plink2("plink2",
      "--pfile", pgen_prefix,
      "--chr",   chr,
      "--recode", "vcf", "bgz",
      "--out", sub("\\.vcf\\.gz$", "", tmp_vcf),
      "--silent"
    )

    args <- c(
      "--input",      tmp_vcf,
      "--map",        gmap,
      "--output",     out_vcf,
      "--region",     glue::glue("chr{chr}"),
      "--thread",     ncores,
      "--pbwt-depth", pbwt_depth
    )
    if (!is.null(scaffold)) {
      args <- c(args, "--scaffold", scaffold,
                "--scaffold-region", glue::glue("chr{chr}"))
    }

    ret <- system2(shapeit4, args, stdout = FALSE, stderr = FALSE)
    if (ret != 0) {
      cli::cli_warn("SHAPEIT4 exited with code {ret} on chr{chr}.")
    } else {
      # Index the output VCF
      system2("tabix", c("-p", "vcf", out_vcf), stdout = FALSE, stderr = FALSE)
      phased_vcfs[as.character(chr)] <- out_vcf
      cli::cli_alert_success("chr{chr} phased → {.path {out_vcf}}")
    }
    unlink(tmp_vcf)
  }

  cli::cli_alert_success("Phasing complete: {length(phased_vcfs)} chromosome(s).")
  invisible(phased_vcfs)
}


#' Infer local ancestry with RFMix2
#'
#' Runs RFMix2 per chromosome on phased study VCFs against a reference panel
#' of labelled ancestry populations. Outputs per-SNP ancestry posterior
#' probability files (\code{.fb.tsv}) which are consumed by
#' \code{\link{gf_load_local_anc}}.
#'
#' @param phased_vcfs    Named character vector of phased study VCF paths,
#'   named by chromosome (output of \code{\link{gf_phase}}).
#' @param ref_vcf        Character. Path to phased reference panel VCF
#'   (e.g., 1KGP3 phased). Must overlap study variants.
#' @param sample_map     Character. Path to 2-column TSV: sample_id, population.
#'   Reference samples only; populations should match the K classes you want
#'   (e.g., AFR, EUR, NAT for Latino cohort).
#' @param outdir         Character. Output directory for RFMix2 output files.
#' @param gmap_dir       Character. Genetic map directory (same as used for
#'   phasing, or Hapi-UR format accepted by RFMix2).
#' @param rfmix2         Character. Path to rfmix2 binary. Default \code{"rfmix2"}.
#' @param em_iterations  Integer. EM iterations. Default 1 (fast; use 5 for
#'   publication-quality runs).
#' @param n_ref_pops     Integer. Expected number of reference populations K.
#'   Used for validation only; inferred from sample_map if NULL.
#' @param chroms         Integer vector. Chromosomes to process. Default same
#'   as names of \code{phased_vcfs}.
#'
#' @return Invisibly returns a named list of output file prefixes per chromosome.
#' @export
#'
#' @examples
#' \dontrun{
#' rfmix_out <- gf_local_ancestry(
#'   phased_vcfs = phased,
#'   ref_vcf     = "data/1kg/kg3_phased.vcf.gz",
#'   sample_map  = "data/1kg/rfmix_sample_map.txt",
#'   outdir      = "data/rfmix",
#'   gmap_dir    = "data/genetic_maps",
#'   chroms      = 21:22
#' )
#' }
gf_local_ancestry <- function(phased_vcfs,
                               ref_vcf,
                               sample_map,
                               outdir        = "data/rfmix",
                               gmap_dir,
                               rfmix2        = "rfmix2",
                               em_iterations = 1L,
                               n_ref_pops    = NULL,
                               chroms        = as.integer(names(phased_vcfs))) {

  .check_binary(rfmix2, "https://github.com/slowkoni/rfmix")
  fs::dir_create(outdir)
  cli::cli_h1("Local Ancestry Inference — RFMix2")

  # Validate reference population count
  if (is.null(n_ref_pops)) {
    smap <- data.table::fread(sample_map, header = FALSE, nThread = 2L)
    n_ref_pops <- length(unique(smap[[2]]))
  }
  cli::cli_alert_info("Reference populations (K): {n_ref_pops}")

  out_prefixes <- list()
  for (chr in chroms) {
    chr_c <- as.character(chr)
    if (!chr_c %in% names(phased_vcfs)) {
      cli::cli_warn("No phased VCF for chr{chr}. Skipping.")
      next
    }
    cli::cli_alert_info("Running RFMix2 on chr{chr}...")
    gmap <- fs::dir_ls(gmap_dir,
                       regexp = glue::glue("chr{chr}[._].*\\.txt(\\.gz)?$"))[1]
    if (is.na(gmap)) {
      cli::cli_warn("No genetic map for chr{chr}. Skipping.")
      next
    }
    out_pfx <- file.path(outdir, glue::glue("chr{chr}"))

    args <- c(
      "-f", phased_vcfs[chr_c],
      "-r", ref_vcf,
      "-m", sample_map,
      "-g", gmap,
      "-o", out_pfx,
      "--chromosome",   glue::glue("chr{chr}"),
      "--n-iterations", em_iterations,
      "--reanalyze-reference"   # ensures reference haplotypes contribute
    )
    ret <- system2(rfmix2, args, stdout = FALSE, stderr = FALSE)
    if (ret != 0) {
      cli::cli_warn("RFMix2 exited with code {ret} on chr{chr}.")
    } else {
      out_prefixes[[chr_c]] <- out_pfx
      cli::cli_alert_success("chr{chr} done → {.path {out_pfx}}.fb.tsv")
    }
  }

  cli::cli_alert_success(
    "Local ancestry inference complete: {length(out_prefixes)} chromosome(s)."
  )
  invisible(out_prefixes)
}


#' Load and align RFMix2 local ancestry posteriors to dosage SNP order
#'
#' Parses RFMix2 \code{.fb.tsv} (forward-backward posterior) files and aligns
#' the per-SNP, per-haplotype ancestry posteriors to the SNP ordering in the
#' dosage matrix. The result is a 4-dimensional array ready for use as the
#' \code{local_anc} input to \code{\link{gf_train}} and
#' \code{\link{gf_predict}}.
#'
#' @details
#' RFMix2 \code{.fb.tsv} files have one row per SNP and columns
#' \code{:::hap0.pop0, :::hap0.pop1, ..., :::hap1.pop0, :::hap1.pop1, ...}
#' giving the posterior probability of each ancestry for each haplotype of
#' each sample. This function:
#' \enumerate{
#'   \item Reads all per-chromosome \code{.fb.tsv} files.
#'   \item Averages the two haplotypes per sample to produce diploid posteriors
#'     (or optionally retains both haplotypes for phased-aware models).
#'   \item Intersects RFMix2 SNPs with dosage matrix SNPs by position.
#'   \item For dosage SNPs absent from RFMix2 output (low-density regions),
#'     fills by nearest-SNP interpolation within each ancestry segment.
#'   \item For samples that are flagged as non-admixed (entropy below
#'     \code{entropy_threshold}), falls back to broadcasting their global
#'     ancestry vector — avoiding noise from uninformative local inference.
#' }
#'
#' @param rfmix_prefixes Named list of RFMix2 output prefixes per chromosome
#'   (output of \code{\link{gf_local_ancestry}}).
#' @param snp_ids        Character vector (M). SNP IDs in dosage matrix order.
#' @param snp_positions  Integer vector (M). Physical positions of dosage SNPs.
#' @param snp_chroms     Integer vector (M). Chromosome of each dosage SNP.
#' @param sample_ids     Character vector (N). Sample IDs in dosage matrix order.
#' @param global_anc     Numeric matrix (N × K). Global ancestry probabilities
#'   from \code{\link{gf_ancestry}}, used as fallback for non-admixed samples.
#' @param keep_haplotypes Logical. If TRUE, return a (N × M × 2 × K) array
#'   retaining haplotype resolution. If FALSE (default), average haplotypes
#'   to return (N × M × K).
#' @param entropy_threshold Numeric. Samples with global ancestry entropy
#'   below this value are treated as non-admixed and filled with global
#'   ancestry vector. Default 0.3 (roughly <10\% minor ancestry component).
#' @param n_cores        Integer. Parallel cores for parsing. Default 2.
#'
#' @return A numeric array of shape \code{(N, M, K)} — or \code{(N, M, 2, K)}
#'   if \code{keep_haplotypes = TRUE} — where N = samples, M = SNPs, K = ref
#'   populations. Dimnames are set on all axes.
#' @export
#'
#' @examples
#' \dontrun{
#' local_anc_arr <- gf_load_local_anc(
#'   rfmix_prefixes = rfmix_out,
#'   snp_ids        = colnames(dosage_matrix),
#'   snp_positions  = snp_pos_vec,
#'   snp_chroms     = snp_chr_vec,
#'   sample_ids     = rownames(dosage_matrix),
#'   global_anc     = anc$anc_probs
#' )
#' dim(local_anc_arr)  # N × M × K
#' }
gf_load_local_anc <- function(rfmix_prefixes,
                               snp_ids,
                               snp_positions,
                               snp_chroms,
                               sample_ids,
                               global_anc,
                               keep_haplotypes     = FALSE,
                               entropy_threshold   = 0.3,
                               n_cores             = 2L) {
  cli::cli_h1("Loading RFMix2 Local Ancestry Posteriors")

  N <- length(sample_ids)
  M <- length(snp_ids)

  # Infer K from first available .fb.tsv header
  first_pfx <- rfmix_prefixes[[1]]
  fb_header  <- data.table::fread(paste0(first_pfx, ".fb.tsv"),
                                   nrows = 0L, nThread = 2L)
  # Column pattern: "#sample:::hap0:::AFR" etc.
  # Extract population names from the header
  hdr       <- names(fb_header)
  pop_names <- unique(sub(".*:::hap[01]:::", "", hdr[grepl(":::", hdr)]))
  K         <- length(pop_names)
  cli::cli_alert_info("K = {K} reference populations: {paste(pop_names, collapse=', ')}")

  # Pre-allocate output array (N × M × K), initialise with global ancestry
  out <- array(NA_real_, dim = c(N, M, K),
               dimnames = list(sample_ids, snp_ids, pop_names))

  # Fill with global ancestry as baseline (broadcast per SNP)
  # Ensures non-admixed / uncovered SNPs have sensible values
  if (!is.null(global_anc)) {
    k_match <- match(pop_names, colnames(global_anc))
    valid_k <- !is.na(k_match)
    for (m in seq_len(M)) {
      out[, m, valid_k] <- global_anc[, k_match[valid_k], drop = FALSE]
    }
    if (any(!valid_k)) {
      cli::cli_warn(
        "Some RFMix2 populations not in global ancestry: {pop_names[!valid_k]}"
      )
    }
  }

  # Identify admixed samples (will receive local ancestry values)
  if (!is.null(global_anc)) {
    entropy <- -.rowSums(
      global_anc * log(global_anc + 1e-9),
      nrow(global_anc), ncol(global_anc)
    )
    admixed_mask <- entropy >= entropy_threshold
    cli::cli_alert_info(
      "{sum(admixed_mask)} admixed samples will use local ancestry; "
      "{sum(!admixed_mask)} non-admixed samples will use global ancestry fallback."
    )
  } else {
    admixed_mask <- rep(TRUE, N)
  }

  # ── Parse per-chromosome .fb.tsv files ────────────────────────────────────
  chr_list <- as.integer(names(rfmix_prefixes))
  pb <- cli::cli_progress_bar("Parsing chromosomes", total = length(chr_list))

  for (chr in chr_list) {
    chr_c  <- as.character(chr)
    fb_file <- paste0(rfmix_prefixes[[chr_c]], ".fb.tsv")
    if (!file.exists(fb_file)) {
      cli::cli_warn("Missing: {.path {fb_file}}. Skipping chr{chr}.")
      cli::cli_progress_update()
      next
    }

    # Which dosage SNPs are on this chromosome?
    chr_mask <- snp_chroms == chr
    if (!any(chr_mask)) {
      cli::cli_progress_update()
      next
    }

    fb <- data.table::fread(fb_file, nThread = n_cores)

    # RFMix2 .fb.tsv columns:
    # #chrom, pos, genetic_pos, n_snps_window,
    # then: SampleID:::hap0:::Pop0, SampleID:::hap0:::Pop1, ...
    #       SampleID:::hap1:::Pop0, SampleID:::hap1:::Pop1, ...
    pos_col <- if ("#pos" %in% names(fb)) "#pos" else "pos"
    rfmix_pos <- fb[[pos_col]]

    # Match RFMix2 positions to dosage SNP positions on this chr
    dosage_pos_chr <- snp_positions[chr_mask]
    dosage_idx_chr <- which(chr_mask)  # global indices in M

    # Nearest-SNP interpolation: for each dosage SNP, find closest RFMix2 SNP
    # This handles SNPs absent from RFMix2 output (low-density / non-genotyped)
    nearest_rfmix <- .nearest_idx(dosage_pos_chr, rfmix_pos)

    # For each sample, extract hap0 + hap1 posteriors and average
    for (n_idx in which(admixed_mask)) {
      sid <- sample_ids[n_idx]
      for (k_idx in seq_len(K)) {
        pop <- pop_names[k_idx]
        col_hap0 <- paste0(sid, ":::hap0:::", pop)
        col_hap1 <- paste0(sid, ":::hap1:::", pop)
        if (!col_hap0 %in% names(fb)) next  # sample not in this RFMix2 run

        h0 <- as.numeric(fb[[col_hap0]])
        h1 <- as.numeric(fb[[col_hap1]])
        # Diploid average of the two haplotypes
        dip <- (h0 + h1) / 2.0
        # Assign to dosage SNP positions via nearest-SNP mapping
        out[n_idx, dosage_idx_chr, k_idx] <- dip[nearest_rfmix]
      }
    }
    cli::cli_progress_update()
  }
  cli::cli_progress_done()

  # Normalise rows to sum to 1 (floating point drift after averaging)
  row_sums <- apply(out, c(1, 2), sum)
  for (k in seq_len(K)) {
    out[, , k] <- out[, , k] / pmax(row_sums, 1e-8)
  }

  cli::cli_alert_success(
    "Local ancestry array ready: [{N} samples \u00d7 {M} SNPs \u00d7 {K} populations]"
  )
  out  # (N × M × K)
}


#' Summary statistics for a local ancestry array
#'
#' Reports mean local ancestry per population per chromosome arm, and flags
#' samples with unusually uniform local ancestry (possible inference failure).
#'
#' @param local_anc_arr Array (N × M × K) from \code{\link{gf_load_local_anc}}.
#' @param snp_chroms    Integer vector (M). Chromosome per SNP.
#' @param pop_labels    Character vector (N). Known superpopulation per sample.
#'
#' @return A data frame with columns pop, chromosome, and mean ancestry per K.
#' @export
gf_local_anc_summary <- function(local_anc_arr, snp_chroms, pop_labels) {
  N   <- dim(local_anc_arr)[1]
  K   <- dim(local_anc_arr)[3]
  pops <- dimnames(local_anc_arr)[[3]]
  chrs <- sort(unique(snp_chroms))

  rows <- list()
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
      rows[[length(rows) + 1]] <- c(
        list(labelled_pop = pop, chr = chr), row
      )
    }
  }
  do.call(rbind, lapply(rows, as.data.frame))
}


# ─── Internal helpers ──────────────────────────────────────────────────────────

.check_binary <- function(bin, url) {
  if (Sys.which(bin) == "" && !file.exists(bin)) {
    cli::cli_abort(c(
      "Binary {.val {bin}} not found on PATH.",
      "i" = "Install from: {url}"
    ))
  }
}

# Nearest-index lookup: for each query position, find index in ref vector
.nearest_idx <- function(query_pos, ref_pos) {
  # Returns integer vector of same length as query_pos
  vapply(query_pos, function(qp) {
    which.min(abs(ref_pos - qp))
  }, integer(1L))
}
