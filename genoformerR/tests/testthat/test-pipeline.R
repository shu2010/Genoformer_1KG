# genoformerR/tests/testthat/test-pipeline.R  (v0.2.0)

# ── Fixtures ───────────────────────────────────────────────────────────────────
make_tiny_dosage <- function(n = 40, m = 50, seed = 1L) {
  set.seed(seed)
  mat <- matrix(sample(0:2, n * m, replace = TRUE, prob = c(.5,.3,.2)),
                nrow = n, ncol = m,
                dimnames = list(
                  paste0("sample_", seq_len(n)),
                  paste0("chr1_", seq_len(m) * 1000, "_A_T")
                ))
  mat
}

make_tiny_panel <- function(n = 40) {
  pops <- c("AFR","AMR","EAS","EUR","SAS")
  data.frame(
    sample    = paste0("sample_", seq_len(n)),
    pop       = sample(pops, n, replace = TRUE),
    super_pop = sample(pops, n, replace = TRUE),
    gender    = sample(c("male","female"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

make_tiny_local_anc <- function(n = 40, m = 50, k = 5, seed = 2L) {
  set.seed(seed)
  # Random Dirichlet-like posteriors: each (n,m) row sums to 1
  raw   <- array(stats::rexp(n * m * k), dim = c(n, m, k))
  norms <- apply(raw, c(1, 2), sum)
  for (ki in seq_len(k)) raw[, , ki] <- raw[, , ki] / norms
  dimnames(raw) <- list(
    paste0("sample_", seq_len(n)),
    paste0("chr1_", seq_len(m) * 1000, "_A_T"),
    c("AFR","AMR","EAS","EUR","SAS")
  )
  raw
}

# ── Ancestry tests (pure R) ────────────────────────────────────────────────────
test_that("gf_ancestry (method='r') returns correct structure", {
  skip_if_not_installed("ranger")
  geno  <- make_tiny_dosage()
  panel <- make_tiny_panel()
  tmp   <- tempfile(fileext = ".panel")
  data.table::fwrite(panel, tmp, sep = "\t")

  anc <- gf_ancestry(dosage_mat = geno, panel_file = tmp,
                     n_pcs = 5, method = "r", n_trees = 50)

  expect_s3_class(anc, "gf_ancestry")
  expect_equal(nrow(anc$pcs), 40)
  expect_equal(ncol(anc$pcs),  5)
  expect_equal(nrow(anc$anc_probs), 40)
  expect_equal(ncol(anc$anc_probs),  5)
  expect_true(all(rowSums(anc$anc_probs) > 0.99))
  expect_equal(length(anc$entropy), 40)
})

test_that("plot.gf_ancestry returns a ggplot", {
  skip_if_not_installed("ranger")
  geno  <- make_tiny_dosage()
  panel <- make_tiny_panel()
  tmp   <- tempfile(fileext = ".panel")
  data.table::fwrite(panel, tmp, sep = "\t")
  anc <- gf_ancestry(geno, tmp, n_pcs = 3, method = "r", n_trees = 20)
  expect_s3_class(plot(anc), "ggplot")
})

# ── Local ancestry array tests ─────────────────────────────────────────────────
test_that("make_tiny_local_anc produces valid posteriors", {
  la <- make_tiny_local_anc()
  expect_equal(dim(la), c(40L, 50L, 5L))
  # Each (n, m) slice should sum to ~1
  row_sums <- apply(la, c(1, 2), sum)
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that("local ancestry array dim names are preserved", {
  la <- make_tiny_local_anc(k = 3)
  expect_equal(dim(la)[3], 3L)
  expect_equal(dimnames(la)[[3]], c("AFR","AMR","EAS"))
})

test_that("gf_local_anc_summary runs without error on synthetic data", {
  la     <- make_tiny_local_anc()
  chroms <- rep(1L, 50)
  labels <- rep(c("AFR","EUR"), each = 20)
  out    <- gf_local_anc_summary(la, chroms, labels)
  expect_s3_class(out, "data.frame")
  expect_true("labelled_pop" %in% names(out))
  expect_true("chr"          %in% names(out))
})

# ── Internal helpers ───────────────────────────────────────────────────────────
test_that(".snp_chroms parses chr:pos format", {
  ids    <- c("chr1_100_A_T", "chr5_200_G_C", "chr22_300_T_A")
  result <- genoformerR:::.snp_chroms(ids)
  expect_equal(result, c(1L, 5L, 22L))
})

test_that(".snp_positions parses chr:pos format", {
  ids <- c("chr1_1000_A_T", "chr5_250000_G_C")
  pos <- genoformerR:::.snp_positions(ids)
  expect_equal(pos, c(1000L, 250000L))
})

test_that(".assign_ld_blocks returns integer vector of correct length", {
  ids    <- paste0("snp", 1:500)
  chroms <- rep(1:5, each = 100)
  blocks <- genoformerR:::.assign_ld_blocks(ids, chroms)
  expect_length(blocks, 500)
  expect_type(blocks, "integer")
})

test_that(".nearest_idx finds correct nearest positions", {
  ref   <- c(100L, 200L, 300L, 500L)
  query <- c(150L, 280L, 450L)
  idx   <- genoformerR:::.nearest_idx(query, ref)
  expect_equal(idx, c(1L, 3L, 4L))
})

# ── Model build (Python, skip if unavailable) ──────────────────────────────────
test_that("gf_build_model with use_local_anc=TRUE builds correctly", {
  skip_if_not(
    tryCatch({ gf_init(); TRUE }, error = function(e) FALSE),
    "Python backend unavailable"
  )
  m <- gf_build_model(n_snps = 100, n_ld_blocks = 20,
                       d_model = 32, n_heads = 4, n_layers = 2,
                       k_pops = 3, d_local = 16, use_local_anc = TRUE)
  expect_s3_class(m, "gf_model")
  expect_true(m$config$use_local_anc)
  expect_equal(m$config$k_pops, 3L)
})

test_that("gf_build_model with use_local_anc=FALSE falls back to global-only", {
  skip_if_not(
    tryCatch({ gf_init(); TRUE }, error = function(e) FALSE),
    "Python backend unavailable"
  )
  m <- gf_build_model(n_snps = 100, n_ld_blocks = 20,
                       d_model = 32, n_heads = 4, n_layers = 2,
                       use_local_anc = FALSE)
  expect_false(m$config$use_local_anc)
  # Global-only model should have fewer parameters
  m2 <- gf_build_model(n_snps = 100, n_ld_blocks = 20,
                        d_model = 32, n_heads = 4, n_layers = 2,
                        k_pops = 5, d_local = 32, use_local_anc = TRUE)
  expect_gt(m2$n_params, m$n_params)
})
