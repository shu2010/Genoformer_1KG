# genoformerR/R/01_data.R
# Step 1–2: 1KGP3 download helpers and PLINK2-based QC wrappers

#' Download 1000 Genomes Phase 3 reference panel
#'
#' Downloads per-chromosome VCFs and the sample panel file from the EBI FTP
#' server. Downloads are parallelised across chromosomes using \code{curl}.
#'
#' @param outdir   Character. Directory for downloaded files.
#' @param chroms   Integer vector. Chromosomes to download. Default 1:22.
#' @param panel_only Logical. If TRUE, only download the sample panel file.
#' @param ncores   Integer. Parallel download workers (uses \code{parallel}).
#'
#' @return Invisibly returns a character vector of downloaded file paths.
#' @export
#'
#' @examples
#' \dontrun{
#' gf_download_1kg(outdir = "data/1kg", chroms = 21:22)
#' }
gf_download_1kg <- function(outdir      = "data/1kg",
                             chroms      = 1:22,
                             panel_only  = FALSE,
                             ncores      = 4L) {
  fs::dir_create(outdir)
  base <- paste0(
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
  )
  panel_url <- paste0(base, "integrated_call_samples_v3.20130502.ALL.panel")

  panel_path <- file.path(outdir, basename(panel_url))
  if (!file.exists(panel_path)) {
    cli::cli_alert_info("Downloading sample panel...")
    download.file(panel_url, panel_path, quiet = TRUE)
  }

  if (panel_only) return(invisible(panel_path))

  vcf_tpl <- paste0(
    base,
    "ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  )
  vcf_urls  <- sprintf(vcf_tpl, chroms)
  vcf_paths <- file.path(outdir, basename(vcf_urls))
  todo      <- vcf_paths[!file.exists(vcf_paths)]

  if (length(todo) > 0) {
    cli::cli_alert_info("Downloading {length(todo)} chromosome VCF(s)...")
    cl <- parallel::makeCluster(min(ncores, length(todo)))
    parallel::parLapply(cl, seq_along(todo), function(i) {
      download.file(vcf_urls[!file.exists(vcf_paths)][i],
                    todo[i], quiet = TRUE, mode = "wb")
    })
    parallel::stopCluster(cl)
  } else {
    cli::cli_alert_info("All VCFs already present.")
  }

  cli::cli_alert_success("1KGP3 data in {.path {outdir}}")
  invisible(c(panel_path, vcf_paths))
}


#' Run PLINK2 QC on 1KGP3 VCFs
#'
#' Converts VCFs to PLINK2 pgen format, merges chromosomes, applies standard
#' GWAS QC filters, performs LD pruning, and intersects variants with a PGS
#' scoring file.
#'
#' @param vcf_dir   Character. Directory containing per-chromosome VCFs.
#' @param pgs_file  Character. Path to PGS Catalog scoring file (.txt.gz).
#' @param outdir    Character. Output directory for pgen files.
#' @param maf       Numeric. Minor allele frequency threshold. Default 0.01.
#' @param hwe       Numeric. HWE p-value threshold. Default 1e-6.
#' @param geno_miss Numeric. SNP missingness threshold. Default 0.02.
#' @param mind_miss Numeric. Sample missingness threshold. Default 0.02.
#' @param ld_window Integer. LD pruning window (variants). Default 200.
#' @param ld_step   Integer. LD pruning step. Default 50.
#' @param ld_r2     Numeric. LD pruning r² threshold. Default 0.3.
#' @param plink2    Character. Path to plink2 binary. Default \code{"plink2"}.
#' @param chroms    Integer vector. Chromosomes to process. Default 1:22.
#'
#' @return A list with elements \code{pgen_prefix}, \code{n_samples},
#'   \code{n_snps_qc}, \code{n_snps_pgs}.
#' @export
gf_qc <- function(vcf_dir,
                   pgs_file,
                   outdir    = "data/plink",
                   maf       = 0.01,
                   hwe       = 1e-6,
                   geno_miss = 0.02,
                   mind_miss = 0.02,
                   ld_window = 200L,
                   ld_step   = 50L,
                   ld_r2     = 0.3,
                   plink2    = "plink2",
                   chroms    = 1:22) {

  fs::dir_create(outdir)
  .check_plink2(plink2)

  cli::cli_h1("Genotype QC Pipeline")

  # ── Step 1: VCF → pgen per chromosome ────────────────────────────────────
  cli::cli_alert_info("Converting VCFs to PLINK2 format...")
  merge_list <- character(0)
  for (chr in chroms) {
    vcf <- fs::dir_ls(vcf_dir, regexp = glue::glue("chr{chr}\\..*\\.vcf\\.gz$"))
    if (length(vcf) == 0) next
    out <- file.path(outdir, glue::glue("chr{chr}_raw"))
    .run_plink2(plink2,
      "--vcf", vcf[1],
      "--make-pgen", "--out", out,
      "--max-alleles", "2", "--rm-dup", "exclude-all",
      "--silent"
    )
    merge_list <- c(merge_list, out)
  }

  # ── Step 2: Merge ─────────────────────────────────────────────────────────
  cli::cli_alert_info("Merging chromosomes...")
  ml_file <- file.path(outdir, "merge_list.txt")
  writeLines(merge_list, ml_file)
  merged  <- file.path(outdir, "kg3_merged")
  .run_plink2(plink2, "--pmerge-list", ml_file, "pfile",
              "--make-pgen", "--out", merged, "--silent")

  # ── Step 3: QC ───────────────────────────────────────────────────────────
  cli::cli_alert_info("Applying QC filters (MAF={maf}, HWE={hwe}, geno={geno_miss})...")
  qc_out <- file.path(outdir, "kg3_qc")
  .run_plink2(plink2,
    "--pfile", merged,
    "--maf",  maf,
    "--hwe",  hwe,
    "--geno", geno_miss,
    "--mind", mind_miss,
    "--make-pgen", "--out", qc_out,
    "--silent"
  )

  # ── Step 4: LD pruning (for ancestry PCA) ─────────────────────────────────
  cli::cli_alert_info("LD pruning for ancestry PCA...")
  prune_out <- file.path(outdir, "pruned_snps")
  .run_plink2(plink2, "--pfile", qc_out,
    "--indep-pairwise", ld_window, ld_step, ld_r2,
    "--out", prune_out, "--silent")
  pruned_out <- file.path(outdir, "kg3_pruned")
  .run_plink2(plink2, "--pfile", qc_out,
    "--extract", paste0(prune_out, ".prune.in"),
    "--make-pgen", "--out", pruned_out, "--silent")

  # ── Step 5: Intersect with PGS scoring file ───────────────────────────────
  cli::cli_alert_info("Intersecting with PGS scoring file...")
  pgs_snps <- .intersect_pgs(pgs_file, paste0(qc_out, ".pvar"))
  snp_list <- file.path(outdir, "pgs_snps.txt")
  writeLines(pgs_snps, snp_list)
  dosage_out <- file.path(outdir, "kg3_pgs_dosage")
  .run_plink2(plink2, "--pfile", qc_out,
    "--extract", snp_list,
    "--export", "A", "--out", dosage_out, "--silent")

  # ── Report ────────────────────────────────────────────────────────────────
  pvar   <- data.table::fread(paste0(qc_out, ".pvar"), nThread = 2L)
  psam   <- data.table::fread(paste0(qc_out, ".psam"), nThread = 2L)
  result <- list(
    pgen_prefix   = qc_out,
    pruned_prefix = pruned_out,
    dosage_raw    = paste0(dosage_out, ".raw"),
    pgs_file      = pgs_file,
    n_samples     = nrow(psam),
    n_snps_qc     = nrow(pvar),
    n_snps_pgs    = length(pgs_snps)
  )
  cli::cli_alert_success(
    "QC complete: {result$n_samples} samples | {result$n_snps_qc} SNPs post-QC | {result$n_snps_pgs} PGS SNPs"
  )
  invisible(result)
}


# ─── Internal helpers ──────────────────────────────────────────────────────────

.check_plink2 <- function(plink2) {
  if (Sys.which(plink2) == "" && !file.exists(plink2)) {
    cli::cli_abort(c(
      "plink2 not found at {.val {plink2}}.",
      "i" = "Install from https://www.cog-genomics.org/plink/2.0/"
    ))
  }
}

.run_plink2 <- function(binary, ...) {
  args <- as.character(c(...))
  ret  <- system2(binary, args, stdout = FALSE, stderr = FALSE)
  if (ret != 0) cli::cli_abort("plink2 exited with code {ret}. Check logs.")
  invisible(ret)
}

.intersect_pgs <- function(pgs_file, pvar_file) {
  pgs  <- data.table::fread(pgs_file, skip = "#", nThread = 2L)
  pvar <- data.table::fread(pvar_file, nThread = 2L)
  # Try rsID column names used by PGS Catalog
  id_col <- intersect(c("rsID", "variant_id", "ID", "snp"), names(pgs))[1]
  if (is.na(id_col)) cli::cli_abort("Cannot find variant ID column in PGS file.")
  intersect(pgs[[id_col]], pvar$ID)
}
