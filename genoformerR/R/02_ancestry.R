# genoformerR/R/02_ancestry.R
# Step 3: PCA-based ancestry inference with soft superpopulation probabilities

#' Infer ancestry from pruned genotypes using PCA + gradient-boosted classifier
#'
#' Computes top principal components from LD-pruned genotype data (using
#' \code{prcomp} natively or the Python sklearn backend), trains a gradient-
#' boosted classifier on known superpopulation labels (AFR/AMR/EAS/EUR/SAS),
#' and returns soft (probabilistic) ancestry proportions for all samples.
#' These proportions serve as conditional inputs to the GenoFormer transformer.
#'
#' @param dosage_mat  Numeric matrix (N × M) of dosage values, OR a path to a
#'   PLINK2 \code{.raw} dosage file. Samples in rows, SNPs in columns.
#' @param panel_file  Character. Path to 1KGP3 sample panel
#'   (\code{integrated_call_samples_v3.20130502.ALL.panel}).
#' @param n_pcs       Integer. Number of PCs to compute. Default 20.
#' @param method      Character. \code{"r"} (native prcomp + ranger) or
#'   \code{"python"} (sklearn via reticulate).
#' @param n_trees     Integer. Number of GBM trees. Default 200.
#' @param seed        Integer. Random seed. Default 42.
#'
#' @return A list with:
#'   \describe{
#'     \item{\code{pcs}}{Numeric matrix (N × n_pcs) of principal components.}
#'     \item{\code{anc_probs}}{Numeric matrix (N × 5) of soft ancestry probabilities
#'       (columns: AFR, AMR, EAS, EUR, SAS).}
#'     \item{\code{predicted_pop}}{Character vector of MAP superpopulation per sample.}
#'     \item{\code{entropy}}{Numeric vector of admixture entropy per sample.}
#'     \item{\code{pca_obj}}{The PCA object (prcomp or sklearn PCA).}
#'     \item{\code{clf_obj}}{The trained classifier object.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' anc <- gf_ancestry(
#'   dosage_mat = "data/plink/kg3_pruned.raw",
#'   panel_file = "data/1kg/integrated_call_samples_v3.20130502.ALL.panel"
#' )
#' plot(anc)
#' }
gf_ancestry <- function(dosage_mat,
                         panel_file,
                         n_pcs    = 20L,
                         method   = c("r", "python"),
                         n_trees  = 200L,
                         seed     = 42L) {
  method <- match.arg(method)
  set.seed(seed)

  # ── Load dosage matrix ─────────────────────────────────────────────────────
  if (is.character(dosage_mat)) {
    cli::cli_alert_info("Reading dosage matrix from {.path {dosage_mat}}...")
    raw <- data.table::fread(dosage_mat, nThread = 4L)
    sample_ids <- raw$IID
    geno <- as.matrix(raw[, -(1:6)])  # drop FID IID PAT MAT SEX PHENOTYPE
    geno[is.na(geno)] <- 0
  } else {
    stopifnot(is.matrix(dosage_mat))
    geno       <- dosage_mat
    sample_ids <- rownames(geno)
  }

  # ── Load panel ─────────────────────────────────────────────────────────────
  panel <- data.table::fread(panel_file, nThread = 2L)
  # Harmonise: 1KGP3 panel columns: sample, pop, super_pop, gender
  if (!"super_pop" %in% names(panel)) {
    cli::cli_abort("Panel file must contain a 'super_pop' column.")
  }
  # Match to dosage sample order
  idx    <- match(sample_ids, panel$sample)
  labels <- panel$super_pop[idx]
  pops   <- c("AFR", "AMR", "EAS", "EUR", "SAS")

  if (method == "r") {
    result <- .ancestry_r(geno, labels, pops, n_pcs, n_trees, seed)
  } else {
    .check_init()
    result <- .ancestry_python(geno, labels, pops, n_pcs, n_trees, seed)
  }

  result$sample_ids <- sample_ids
  result$labels     <- labels
  class(result)     <- c("gf_ancestry", "list")

  cli::cli_alert_success(
    "Ancestry inference complete. Mean admixture entropy: {.val {round(mean(result$entropy), 3)}}"
  )
  result
}


#' @export
print.gf_ancestry <- function(x, ...) {
  pops <- c("AFR", "AMR", "EAS", "EUR", "SAS")
  cli::cli_h1("GenoFormer Ancestry Inference")
  cli::cli_alert_info("{length(x$sample_ids)} samples | {ncol(x$pcs)} PCs")
  cat("\nMean soft ancestry proportions by labelled superpopulation:\n")
  for (p in pops) {
    mask <- x$labels == p
    if (!any(mask, na.rm = TRUE)) next
    mean_probs <- round(colMeans(x$anc_probs[mask, , drop = FALSE]), 3)
    cat(sprintf("  %-5s  %s\n", p,
        paste(sprintf("%s=%.3f", pops, mean_probs), collapse = "  ")))
  }
  cat(sprintf("\nMean admixture entropy: %.3f\n", mean(x$entropy)))
  invisible(x)
}


#' Plot ancestry PCA coloured by superpopulation
#'
#' @param x   A \code{gf_ancestry} object.
#' @param pcs Integer vector of length 2. Which PCs to plot. Default \code{c(1,2)}.
#' @param ...  Unused.
#'
#' @return A \code{ggplot2} object.
#' @export
plot.gf_ancestry <- function(x, pcs = c(1L, 2L), ...) {
  pal <- c(AFR = "#ff6b6b", AMR = "#ffd700", EAS = "#69ff47",
           EUR = "#00e5ff", SAS = "#c77dff")
  df <- data.frame(
    PC_x = x$pcs[, pcs[1]],
    PC_y = x$pcs[, pcs[2]],
    pop  = factor(x$labels, levels = names(pal))
  )
  ggplot2::ggplot(df, ggplot2::aes(PC_x, PC_y, colour = pop)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.8) +
    ggplot2::scale_colour_manual(values = pal, na.value = "grey60") +
    ggplot2::labs(
      title   = "1KGP3 Ancestry PCA",
      x       = paste0("PC", pcs[1]),
      y       = paste0("PC", pcs[2]),
      colour  = "Superpopulation"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = "#0a0a14", colour = NA),
      panel.background = ggplot2::element_rect(fill = "#0f0f1e", colour = NA),
      text             = ggplot2::element_text(colour = "#e0e0ff"),
      axis.text        = ggplot2::element_text(colour = "#8888bb"),
      panel.grid       = ggplot2::element_line(colour = "#1a1a3a"),
      legend.background = ggplot2::element_rect(fill = "#0a0a14", colour = NA)
    )
}


# ─── Internal implementations ──────────────────────────────────────────────────

.ancestry_r <- function(geno, labels, pops, n_pcs, n_trees, seed) {
  if (!requireNamespace("ranger", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg ranger} required for method='r'. Install with: install.packages('ranger')")
  }
  cli::cli_alert_info("Running PCA (prcomp)...")
  pca <- stats::prcomp(geno, center = TRUE, scale. = FALSE, rank. = n_pcs)
  pcs <- pca$x[, seq_len(n_pcs)]

  cli::cli_alert_info("Training gradient-boosted ancestry classifier (ranger)...")
  df_train <- as.data.frame(pcs)
  df_train$pop <- factor(labels, levels = pops)
  # Drop samples with missing labels
  df_train <- df_train[!is.na(df_train$pop), ]

  clf <- ranger::ranger(
    pop ~ .,
    data           = df_train,
    num.trees      = n_trees,
    probability    = TRUE,
    seed           = seed,
    num.threads    = parallel::detectCores() - 1L
  )
  anc_probs <- stats::predict(clf, data = as.data.frame(pcs))$predictions
  colnames(anc_probs) <- levels(df_train$pop)

  predicted_pop <- pops[max.col(anc_probs)]
  entropy       <- -.rowSums(anc_probs * log(anc_probs + 1e-9), nrow(anc_probs), ncol(anc_probs))

  list(pcs = pcs, anc_probs = anc_probs,
       predicted_pop = predicted_pop, entropy = entropy,
       pca_obj = pca, clf_obj = clf)
}


.ancestry_python <- function(geno, labels, pops, n_pcs, n_trees, seed) {
  .check_init()
  cli::cli_alert_info("Running PCA (sklearn)...")
  sk_dec <- reticulate::import("sklearn.decomposition", convert = FALSE)
  sk_ens <- reticulate::import("sklearn.ensemble",      convert = FALSE)
  sk_pre <- reticulate::import("sklearn.preprocessing", convert = FALSE)
  np     <- .gf_env$np

  geno_py  <- np$array(geno,   dtype = "float32")
  label_py <- np$array(labels, dtype = "str")

  pca     <- sk_dec$PCA(n_components = as.integer(n_pcs), random_state = as.integer(seed))
  pcs_py  <- pca$fit_transform(geno_py)

  le      <- sk_pre$LabelEncoder()
  y       <- le$fit_transform(label_py)

  clf     <- sk_ens$GradientBoostingClassifier(
    n_estimators = as.integer(n_trees),
    max_depth    = 4L,
    random_state = as.integer(seed)
  )
  clf$fit(pcs_py, y)
  anc_probs <- reticulate::py_to_r(clf$predict_proba(pcs_py))

  pcs_r         <- reticulate::py_to_r(pcs_py)
  colnames(anc_probs) <- reticulate::py_to_r(le$classes_)
  predicted_pop <- colnames(anc_probs)[max.col(anc_probs)]
  entropy       <- -.rowSums(anc_probs * log(anc_probs + 1e-9), nrow(anc_probs), ncol(anc_probs))

  list(pcs = pcs_r, anc_probs = anc_probs,
       predicted_pop = predicted_pop, entropy = entropy,
       pca_obj = pca, clf_obj = clf)
}
