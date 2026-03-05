# genoformerR/R/zzz.R
# Package initialisation — Python environment bootstrap via reticulate

.gf_env <- new.env(parent = emptyenv())

#' @importFrom reticulate py_module_available import use_condaenv conda_create
#'   conda_install py_run_file
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_danger cli_h1
NULL

.onLoad <- function(libname, pkgname) {
  # Defer Python import until first use; store module handles in .gf_env
  .gf_env$initialized <- FALSE
}

#' Initialise the GenoFormer Python backend
#'
#' Creates (or reuses) a conda environment named \code{"genoformer"} and
#' imports all required Python modules. Call this once at the start of a
#' session before any modelling functions.
#'
#' @param envname Character. Conda environment name. Default \code{"genoformer"}.
#' @param python  Character or NULL. Path to a specific Python binary; if NULL,
#'   reticulate selects automatically.
#' @param reinstall Logical. Force reinstallation of Python dependencies.
#'
#' @return Invisibly returns a named list of imported Python module handles.
#' @export
#'
#' @examples
#' \dontrun{
#' gf_init()
#' }
gf_init <- function(envname = "genoformer", python = NULL, reinstall = FALSE) {
  cli::cli_h1("GenoFormer Initialisation")

  # ── 1. Conda environment ──────────────────────────────────────────────────
  envs <- tryCatch(reticulate::conda_list()$name, error = function(e) character(0))
  if (!(envname %in% envs) || reinstall) {
    cli::cli_alert_info("Creating conda environment '{envname}'...")
    reticulate::conda_create(envname, python_version = "3.10")
    reticulate::conda_install(
      envname,
      packages = c(
        "pytorch", "torchvision", "torchaudio",  # from pytorch channel
        "numpy", "pandas", "scikit-learn",
        "scipy", "matplotlib", "tqdm"
      ),
      channel  = c("pytorch", "conda-forge"),
      pip      = FALSE
    )
    # Install plink-style IO helpers via pip
    reticulate::py_install(
      c("pyreadr", "h5py"),
      envname = envname, pip = TRUE
    )
    cli::cli_alert_success("Environment '{envname}' ready.")
  }

  if (!is.null(python)) {
    reticulate::use_python(python, required = TRUE)
  } else {
    reticulate::use_condaenv(envname, required = TRUE)
  }

  # ── 2. Import Python modules ──────────────────────────────────────────────
  cli::cli_alert_info("Importing Python modules...")
  .gf_env$np      <- reticulate::import("numpy",        convert = FALSE)
  .gf_env$pd      <- reticulate::import("pandas",       convert = FALSE)
  .gf_env$torch   <- reticulate::import("torch",        convert = FALSE)
  .gf_env$nn      <- reticulate::import("torch.nn",     convert = FALSE)
  .gf_env$sklearn <- reticulate::import("sklearn",      convert = FALSE)
  .gf_env$scipy   <- reticulate::import("scipy",        convert = FALSE)

  # ── 3. Source Python pipeline scripts bundled with the package ─────────────
  pkg_py <- system.file("python", package = "genoformerR")
  for (script in c("model.py", "train.py", "evaluate.py")) {
    path <- file.path(pkg_py, script)
    if (file.exists(path)) reticulate::py_run_file(path)
  }

  .gf_env$initialized <- TRUE
  cli::cli_alert_success("GenoFormer backend ready (PyTorch {.val {reticulate::py_to_r(.gf_env$torch$`__version__`)}})")
  invisible(.gf_env)
}

# Internal: assert backend is ready
.check_init <- function() {
  if (!isTRUE(.gf_env$initialized)) {
    cli::cli_abort("Run {.code gf_init()} first to initialise the Python backend.")
  }
}
