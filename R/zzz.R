################################################################################
# Package Loading and Initialization
################################################################################

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Set TBB environment variable (from bmm_core.R)
  if (Sys.getenv("TBB_CXX_TYPE") == "") {
    Sys.setenv(TBB_CXX_TYPE = "gcc")
  }
  
  # Set parallel cores (from bmm_core.R)
  options(mc.cores = parallel::detectCores())
  # Check TwoSampleMR
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }
    remotes::install_github("MRCIEU/TwoSampleMR")
  }
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "\n",
    "╔════════════════════════════════════════════════════════╗\n",
    "║  BMM: Bayesian Mixture Model for Mendelian            ║\n",
    "║       Randomization                                    ║\n",
    "║                                                        ║\n",
    "║  Version 1.0.0                                         ║\n",
    "╚════════════════════════════════════════════════════════╝\n",
    "\n",
    "⚠️  IMPORTANT - Environment setup:\n",
    "   Sys.setenv(TBB_CXX_TYPE = \"gcc\")\n",
    "   (Already set for this session)\n",
    "\n",
    "Add to ~/.Rprofile for permanent setup:\n",
    "  if (interactive()) Sys.setenv(TBB_CXX_TYPE = \"gcc\")\n",
    "\n",
    "Quick Start:\n",
    "  ?run_mr_model_comparison\n",
    "  ?harmonize_for_mr\n",
    "  install_bmm_dependencies()  # Install cmdstanr + CmdStan\n",
    "\n"
  )
  
  # Check cmdstanr and CmdStan installation
  check_stan_setup()
}

#' Check Stan Setup
#' @keywords internal
check_stan_setup <- function() {
  # Check cmdstanr package
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    packageStartupMessage(
      "⚠️  cmdstanr package not installed\n\n",
      "Install with:\n",
      "  install.packages('cmdstanr',\n",
      "    repos = c('https://stan-dev.r-universe.dev',\n",
      "              getOption('repos')))\n\n",
      "Or use helper:\n",
      "  install_bmm_dependencies()\n"
    )
    return(invisible(FALSE))
  }
  
  # Check CmdStan binary - CORRECT METHOD
  cmdstan_available <- tryCatch({
    version <- cmdstanr::cmdstan_version()
    TRUE
  }, error = function(e) {
    FALSE
  })
  
  if (!cmdstan_available) {
    packageStartupMessage(
      "⚠️  CmdStan not installed\n\n",
      "After installing cmdstanr, run:\n",
      "  cmdstanr::install_cmdstan()\n\n",
      "Or use helper:\n",
      "  install_bmm_dependencies()\n"
    )
    return(invisible(FALSE))
  }
  
  # All OK - get version safely
  version <- tryCatch({
    cmdstanr::cmdstan_version()
  }, error = function(e) {
    "unknown"
  })
  
  packageStartupMessage(
    "✓ cmdstanr and CmdStan are ready\n",
    sprintf("  CmdStan version: %s\n", version)
  )
  
  return(invisible(TRUE))
}