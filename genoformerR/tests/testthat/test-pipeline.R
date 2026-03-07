# genoformerR/tests/testthat/test-pipeline.R  (v0.3.0)

# ── Fixtures ───────────────────────────────────────────────────────────────────
make_tiny_dosage <- function(n = 60, m = 120, seed = 1L) {
  set.seed(seed)
  mat <- matrix(
    sample(0:2, n * m, replace = TRUE, prob = c(.5, .3, .2)),
    nrow = n, ncol = m,
    dimnames = list(
      paste0("sample_", seq_len(n)),
      paste0("chr", rep(1:3, each = m / 3), "_",
             rep(seq(1000, by = 1000, length.out = m / 3), 3), "_A_T")
    )
  )
  mat
}

make_tiny_panel <- function(n = 60) {
  pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  data.frame(
    sample    = paste0("sample_", seq_len(n)),
    pop       = sample(pops, n, replace = TRUE),
    super_pop = sample(pops, n, replace = TRUE),
    gender    = sample(c("male", "female"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

make_tiny_local_anc <- function(n = 60, m = 120, k = 5, seed = 2L) {
  set.seed(seed)
  raw   <- array(stats::rexp(n * m * k), dim = c(n, m, k))
  norms <- apply(raw, c(1, 2), sum)
  for (ki in seq_len(k)) raw[, , ki] <- raw[, , ki] / norms
  dimnames(raw) <- list(
    paste0("sample_", seq_len(n)),
    paste0("chr", rep(1:3, each = m / 3), "_",
           rep(seq(1000, by = 1000, length.out = m / 3), 3), "_A_T"),
    c("AFR", "AMR", "EAS", "EUR", "SAS")
  )
  raw
}

# ── Global ancestry ────────────────────────────────────────────────────────────
test_that("gf_ancestry (method='r') returns correct structure", {
  skip_if_not_installed("ranger")
  geno  <- make_tiny_dosage()
  panel <- make_tiny_panel()
  tmp   <- tempfile(fileext = ".panel")
  data.table::fwrite(panel, tmp, sep = "\t")

  anc <- gf_ancestry(dosage_mat = geno, panel_file = tmp,
                     n_pcs = 5, method = "r", n_trees = 50)
  expect_s3_class(anc, "gf_ancestry")
  expect_equal(nrow(anc$pcs), 60)
  expect_equal(ncol(anc$pcs),  5)
  expect_true(all(rowSums(anc$anc_probs) > 0.99))
})

test_that("plot.gf_ancestry returns ggplot", {
  skip_if_not_installed("ranger")
  geno  <- make_tiny_dosage()
  panel <- make_tiny_panel()
  tmp   <- tempfile(fileext = ".panel"); data.table::fwrite(panel, tmp, sep = "\t")
  anc   <- gf_ancestry(geno, tmp, n_pcs = 3, method = "r", n_trees = 20)
  expect_s3_class(plot(anc), "ggplot")
})

# ── Windowed PCA ───────────────────────────────────────────────────────────────
test_that("gf_windowed_pca returns correct dimensions", {
  geno   <- make_tiny_dosage(n = 30, m = 60)
  chrs   <- .snp_chroms(colnames(geno))
  la     <- gf_windowed_pca(geno, snp_chroms = chrs, k_pops = 3,
                              window_size = 20)
  expect_equal(dim(la), c(30L, 60L, 3L))
})

test_that("gf_windowed_pca rows sum to 1", {
  geno   <- make_tiny_dosage(n = 20, m = 40)
  chrs   <- rep(1L, 40)
  la     <- gf_windowed_pca(geno, snp_chroms = chrs, k_pops = 4,
                              window_size = 10)
  row_sums <- apply(la, c(1, 2), sum)
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that("gf_windowed_pca handles multi-chromosome input", {
  geno   <- make_tiny_dosage()
  chrs   <- .snp_chroms(colnames(geno))
  expect_equal(length(unique(chrs)), 3L)  # fixture has 3 chromosomes
  la     <- gf_windowed_pca(geno, snp_chroms = chrs, k_pops = 5,
                              window_size = 30)
  expect_equal(dim(la)[1:2], c(60L, 120L))
})

test_that("gf_windowed_pca respects k_pops capping at N-1", {
  geno <- make_tiny_dosage(n = 4, m = 20)
  chrs <- rep(1L, 20)
  # k_pops = 10 should be capped to 3 (N-1 = 3)
  la <- gf_windowed_pca(geno, snp_chroms = chrs, k_pops = 10,
                         window_size = 10)
  expect_equal(dim(la)[3], 3L)
})

test_that("gf_windowed_pca with overlapping windows (step < window)", {
  geno <- make_tiny_dosage(n = 20, m = 60)
  chrs <- rep(1L, 60)
  la   <- gf_windowed_pca(geno, snp_chroms = chrs, k_pops = 3,
                           window_size = 20, step_size = 10)
  expect_equal(dim(la), c(20L, 60L, 3L))
  row_sums <- apply(la, c(1, 2), sum)
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

# ── gf_local_ancestry dispatcher ──────────────────────────────────────────────
test_that("gf_local_ancestry(method='windowed_pca') works end-to-end", {
  geno <- make_tiny_dosage(n = 20, m = 40)
  chrs <- rep(1L, 40)
  la   <- gf_local_ancestry(
    method      = "windowed_pca",
    dosage_mat  = geno,
    snp_chroms  = chrs,
    k_pops      = 3L,
    window_size = 15
  )
  expect_equal(dim(la)[3], 3L)
})

test_that("gf_local_ancestry(method='flare') errors without required args", {
  expect_error(
    gf_local_ancestry(method = "flare", dosage_mat = matrix(0, 5, 10)),
    "requires"
  )
})

test_that("gf_local_ancestry(method='rfmix2') errors without required args", {
  expect_error(
    gf_local_ancestry(method = "rfmix2", dosage_mat = matrix(0, 5, 10)),
    "requires"
  )
})

# ── Local ancestry array utilities ────────────────────────────────────────────
test_that("gf_local_anc_summary returns correct structure", {
  la     <- make_tiny_local_anc()
  chrs   <- .snp_chroms(dimnames(la)[[2]])
  labels <- rep(c("AFR", "EUR"), each = 30)
  out    <- gf_local_anc_summary(la, chrs, labels)
  expect_s3_class(out, "data.frame")
  expect_true("labelled_pop" %in% names(out))
  expect_true("chr"          %in% names(out))
})

# ── Internal helpers ───────────────────────────────────────────────────────────
test_that(".snp_chroms parses chr:pos format", {
  ids <- c("chr1_1000_A_T", "chr5_2000_G_C", "chr22_3000_T_A")
  expect_equal(genoformerR:::.snp_chroms(ids), c(1L, 5L, 22L))
})

test_that(".snp_positions parses chr:pos format", {
  ids <- c("chr1_1000_A_T", "chr5_250000_G_C")
  expect_equal(genoformerR:::.snp_positions(ids), c(1000L, 250000L))
})

test_that(".assign_ld_blocks returns integer vector", {
  ids    <- paste0("snp", 1:500)
  chroms <- rep(1:5, each = 100)
  blocks <- genoformerR:::.assign_ld_blocks(ids, chroms)
  expect_length(blocks, 500)
  expect_type(blocks, "integer")
})

test_that(".nearest_idx finds correct nearest positions", {
  ref   <- c(100L, 200L, 300L, 500L)
  query <- c(150L, 280L, 450L)
  expect_equal(genoformerR:::.nearest_idx(query, ref), c(1L, 3L, 4L))
})

# ── Model build (Python, skip if unavailable) ──────────────────────────────────
test_that("gf_build_model with use_local_anc=TRUE builds correctly", {
  skip_if_not(
    tryCatch({ gf_init(); TRUE }, error = function(e) FALSE),
    "Python backend unavailable"
  )
  m <- gf_build_model(n_snps = 100, n_ld_blocks = 20, d_model = 32,
                       n_heads = 4, n_layers = 2, k_pops = 3,
                       d_local = 16, use_local_anc = TRUE)
  expect_s3_class(m, "gf_model")
  expect_true(m$config$use_local_anc)
  expect_equal(m$config$k_pops, 3L)
})

test_that("gf_build_model global-only has fewer params than dual-level", {
  skip_if_not(
    tryCatch({ gf_init(); TRUE }, error = function(e) FALSE),
    "Python backend unavailable"
  )
  m_g  <- gf_build_model(n_snps = 100, n_ld_blocks = 20, d_model = 32,
                           n_heads = 4, n_layers = 2, use_local_anc = FALSE)
  m_dl <- gf_build_model(n_snps = 100, n_ld_blocks = 20, d_model = 32,
                           n_heads = 4, n_layers = 2, k_pops = 5,
                           d_local = 32, use_local_anc = TRUE)
  expect_gt(m_dl$n_params, m_g$n_params)
})

# ── Pipeline argument validation ───────────────────────────────────────────────
test_that("gf_run_pipeline errors when FLARE args missing", {
  skip_if_not(
    tryCatch({ gf_init(); TRUE }, error = function(e) FALSE),
    "Python backend unavailable"
  )
  expect_error(
    gf_run_pipeline(
      dosage_raw       = tempfile(),
      pgs_file         = tempfile(),
      panel_file       = tempfile(),
      phenotype        = numeric(0),
      local_anc_method = "flare"
      # vcf_file, ref_vcf etc. deliberately omitted
    ),
    "requires"
  )
})
