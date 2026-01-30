################################################################################
# BMM Wrapper Functions with roxygen2 documentation
################################################################################

#' Run BMM Analysis
#' 
#' Wrapper for run_mr_model_comparison() from bmm_core.R
#' Performs complete Bayesian mixture model analysis with 3/4-component
#' model comparison and automatic model selection.
#' 
#' @param beta_X Numeric vector of SNP-exposure associations
#' @param beta_Y Numeric vector of SNP-outcome associations  
#' @param se_X Numeric vector of standard errors for beta_X
#' @param se_Y Numeric vector of standard errors for beta_Y
#' @param rho_e Numeric scalar, sample overlap correlation (default: 0)
#' @param compiled_models List with compiled Stan models from compile_stan_models()
#' @param n_chains Integer, number of MCMC chains (default: 4)
#' @param iter Integer, total iterations per chain (default: 1000)
#' @param warmup Integer, warmup iterations (default: 200)
#' @param use_weighted_mode Logical, use Weighted Mode initialization (default: TRUE)
#' @param use_parallel Logical, run models in parallel (default: TRUE)
#' @param ... Additional arguments passed to run_mr_model_comparison()
#' 
#' @return Object of class 'bmm_result' with analysis results
#' 
#' @details
#' This function performs:
#' \enumerate{
#'   \item Fits 3-component model (Valid, Valid+UHP, NULL/UHP)
#'   \item Fits 4-component model (Valid, Valid+UHP, Valid+UHP+CHP, NULL)
#'   \item Compares models via LOO-CV
#'   \item Selects best model
#'   \item Returns final inference
#' }
#' 
#' IMPORTANT: Set Sys.setenv(TBB_CXX_TYPE = "gcc") before running.
#' 
#' @examples
#' \dontrun{
#' # Set environment
#' Sys.setenv(TBB_CXX_TYPE = "gcc")
#' 
#' # Compile models
#' compiled <- compile_stan_models()
#' 
#' # Run analysis
#' results <- run_bmm_analysis(
#'   beta_X = harmonized_data$beta_X,
#'   beta_Y = harmonized_data$beta_Y,
#'   se_X = harmonized_data$se_X,
#'   se_Y = harmonized_data$se_Y,
#'   rho_e = 0,
#'   compiled_models = compiled
#' )
#' 
#' # View results
#' summary(results)
#' }
#' 
#' @seealso \code{\link{compile_stan_models}}, \code{\link{harmonize_for_mr}}
#' 
#' @export
run_bmm_analysis <- function(beta_X, beta_Y, se_X, se_Y, 
                             rho_e = 0,
                             compiled_models = NULL,
                             n_chains = 1,
                             iter = 1000,
                             warmup = 200,
                             use_weighted_mode = TRUE,
                             use_parallel = FALSE,
                             ...) {
  
  # Ensure TBB is set
  if (Sys.getenv("TBB_CXX_TYPE") == "") {
    warning("TBB_CXX_TYPE not set. Setting to 'gcc'...\n",
            "Add Sys.setenv(TBB_CXX_TYPE = 'gcc') to your script.",
            call. = FALSE)
    Sys.setenv(TBB_CXX_TYPE = "gcc")
  }
  
  # Call original function from bmm_core.R
  results <- run_mr_model_comparison(
    beta_X = beta_X,
    beta_Y = beta_Y,
    se_X = se_X,
    se_Y = se_Y,
    rho_e = rho_e,
    compiled_models = compiled_models,
    n_chains = n_chains,
    iter = iter,
    warmup = warmup,
    use_weighted_mode = use_weighted_mode,
    use_parallel = use_parallel,
    ...
  )
  return(results)
}

#' Compile BMM Stan Models
#' 
#' Wrapper for compile_stan_models() with automatic TBB environment setup.
#' Stan models are compiled from source and cannot be pre-compiled.
#' 
#' @param force_recompile Logical, force recompilation even if models exist (default: FALSE)
#' @param set_tbb_env Logical, automatically set TBB_CXX_TYPE environment variable (default: TRUE)
#' 
#' @return List with two compiled cmdstan_model objects:
#' \itemize{
#'   \item model_3comp - 3-component model
#'   \item model_4comp - 4-component model
#' }
#' 
#' @details
#' This function compiles Stan models from source code strings defined in bmm_core.R.
#' Compilation typically takes 1-2 minutes and must be done once per R session.
#' 
#' The models cannot be pre-compiled because Stan models compile to system-specific
#' C++ code. Each user must compile fresh models on their system.
#' 
#' @examples
#' \dontrun{
#' # Compile models (required before analysis)
#' Sys.setenv(TBB_CXX_TYPE = "gcc")
#' compiled <- compile_stan_models()
#' 
#' # Force recompilation
#' compiled <- compile_stan_models(force_recompile = TRUE)
#' }
#' 
#' @seealso \code{\link{run_bmm_analysis}}
#' 
#' @export
compile_bmm_models <- function(force_recompile = FALSE, set_tbb_env = TRUE) {
  
  cat("\n╔════════════════════════════════════════╗\n")
  cat("║  Compiling BMM Stan Models            ║\n")
  cat("╚════════════════════════════════════════╝\n\n")
  
  # Set TBB environment
  if (set_tbb_env) {
    current_tbb <- Sys.getenv("TBB_CXX_TYPE")
    if (current_tbb == "") {
      cat("→ Setting TBB_CXX_TYPE = 'gcc'\n\n")
      Sys.setenv(TBB_CXX_TYPE = "gcc")
    } else {
      cat(sprintf("→ TBB_CXX_TYPE = '%s'\n\n", current_tbb))
    }
  }
  
  cat("⚠️  IMPORTANT:\n")
  cat("  • Models compile from source (takes 1-2 min)\n")
  cat("  • Cannot be pre-compiled (system-specific)\n")
  cat("  • Required once per R session\n\n")
  
  # Call original function from bmm_core.R
  compiled <- compile_stan_models(force_recompile = force_recompile)
  
  cat("\n✓ Models compiled successfully!\n")
  cat("  You can reuse these for multiple analyses.\n\n")
  
  return(compiled)
}

#' Setup BMM Environment
#' 
#' Check and setup BMM requirements including TBB, cmdstanr, and CmdStan.
#' 
#' @return Invisibly returns TRUE if all requirements met, FALSE otherwise
#' 
#' @details
#' Checks:
#' \itemize{
#'   \item TBB_CXX_TYPE environment variable
#'   \item cmdstanr package installation
#'   \item CmdStan installation
#'   \item plink availability (optional, for harmonization)
#' }
#' 
#' @examples
#' \dontrun{
#' # Check environment
#' setup_bmm_environment()
#' }
#' 
#' @export
setup_bmm_environment <- function() {
  cat("\n═══ BMM Environment Setup ═══\n\n")
  
  all_ok <- TRUE
  
  # TBB
  tbb <- Sys.getenv("TBB_CXX_TYPE")
  if (tbb == "") {
    Sys.setenv(TBB_CXX_TYPE = "gcc")
    cat("✓ TBB_CXX_TYPE set to 'gcc'\n")
  } else {
    cat(sprintf("✓ TBB_CXX_TYPE = '%s'\n", tbb))
  }
  
  # cmdstanr
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    cat("✗ cmdstanr not installed\n")
    cat("  Install: install.packages('cmdstanr',\n")
    cat("             repos = c('https://mc-stan.org/r-packages/',\n")
    cat("                       getOption('repos')))\n")
    all_ok <- FALSE
  } else {
    cat("✓ cmdstanr installed\n")
    
    cmdstan_ok <- tryCatch({
      cmdstanr::cmdstan_version()
      TRUE
    }, error = function(e) {
      FALSE
    })
    
    if (!cmdstan_ok) {
      cat("✗ CmdStan not installed\n")
      cat("  Install: cmdstanr::install_cmdstan()\n")
      all_ok <- FALSE
    } else {
      cat("✓ CmdStan installed\n")
      version <- tryCatch({
        cmdstanr::cmdstan_version()
      }, error = function(e) {
        "unknown"
      })
      cat(sprintf("  Version: %s\n", version))
    }
  }
  
  # Plink check (optional)
  plink_paths <- c("/usr/bin/plink", "/usr/local/bin/plink", 
                   "/opt/homebrew/bin/plink")
  plink_found <- FALSE
  for (p in plink_paths) {
    if (file.exists(p)) {
      cat(sprintf("✓ Plink found: %s\n", p))
      plink_found <- TRUE
      break
    }
  }
  if (!plink_found) {
    cat("ℹ Plink not found in common locations\n")
    cat("  (Needed for harmonization with proxy SNP search)\n")
  }
  
  cat("\n")
  if (all_ok) {
    cat("✓ Environment ready!\n\n")
  } else {
    cat("⚠️  Some requirements missing. See above.\n\n")
  }
  
  invisible(all_ok)
}