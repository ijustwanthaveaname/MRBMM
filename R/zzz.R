################################################################################
# Package Loading and Initialization (SAFE VERSION)
################################################################################

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Set TBB environment variable if not already set
  if (Sys.getenv("TBB_CXX_TYPE") == "") {
    Sys.setenv(TBB_CXX_TYPE = "gcc")
  }
  invisible()
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "BMM: Bayesian Mixture Model for Mendelian Randomization\n",
    "Version 0.0.1\n",
    "\n",
    "This package requires external dependencies for full functionality.\n",
    "Run the following after loading the package:\n",
    "\n",
    "  install_bmm_dependencies()   # Install non-CRAN dependencies\n",
    "  check_stan_setup()           # Verify cmdstanr and CmdStan\n",
    "\n",
    "See github.com/ijustwanthaveaname/MRBMM/README.md for usage.\n"
  )
}

################################################################################
# User-invoked helper functions
################################################################################

#' Check Stan setup for BMM
#'
#' Verify that \pkg{cmdstanr} and CmdStan are correctly installed
#' and available for model fitting.
#'
#' This function is intentionally NOT called during package loading.
#'
#' @return Invisibly returns TRUE if CmdStan is available, FALSE otherwise.
#' @export
check_stan_setup <- function() {

  # Check cmdstanr package
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    message(
      "cmdstanr is not installed.\n\n",
      "Install with:\n",
      "  install.packages(\n",
      "    'cmdstanr',\n",
      "    repos = c('https://stan-dev.r-universe.dev', getOption('repos'))\n",
      "  )\n"
    )
    return(invisible(FALSE))
  }

  # Check CmdStan binary availability
  cmdstan_available <- tryCatch({
    cmdstanr::cmdstan_version()
    TRUE
  }, error = function(e) {
    FALSE
  })

  if (!cmdstan_available) {
    message(
      "CmdStan is not installed.\n\n",
      "After installing cmdstanr, run:\n",
      "  cmdstanr::install_cmdstan()\n"
    )
    return(invisible(FALSE))
  }

  version <- tryCatch(
    cmdstanr::cmdstan_version(),
    error = function(e) "unknown"
  )

  message(
    "✓ cmdstanr and CmdStan are ready\n",
    sprintf("  CmdStan version: %s\n", version)
  )

  invisible(TRUE)
}

#' Install non-CRAN dependencies for BMM
#'
#' Install GitHub-hosted and system-level dependencies required
#' for running the full BMM pipeline.
#'
#' This function must be called explicitly by the user.
#'
#' @export
install_bmm_dependencies <- function() {

  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
  }

  github_pkgs <- c(
    "MRCIEU/TwoSampleMR",
    "MRCIEU/ieugwasr"
  )

  for (repo in github_pkgs) {
    pkg <- basename(repo)
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing ", pkg, " from GitHub...")
      remotes::install_github(repo)
    }
  }

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    message("Installing cmdstanr...")
    remotes::install_github("stan-dev/cmdstanr")
  }

  message("✓ BMM dependencies installation complete")
  invisible(TRUE)
}
