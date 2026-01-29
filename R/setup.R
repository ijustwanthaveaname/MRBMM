################################################################################
# BMM Setup and Dependency Installation
################################################################################

#' Check if CmdStan is Installed
#' 
#' @return Logical, TRUE if installed
#' @keywords internal
check_cmdstan_installed <- function() {
  tryCatch({
    cmdstanr::cmdstan_version()
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Install BMM Dependencies
#' 
#' Helper function to install cmdstanr and CmdStan.
#' This is recommended for first-time users.
#' 
#' @param install_cmdstan Logical, also install CmdStan binary (default: TRUE)
#' @param ask Logical, ask for confirmation before installing (default: TRUE)
#' @param cores Integer, number of cores for CmdStan compilation (default: auto-detect)
#' @param verbose Logical, print progress messages (default: TRUE)
#' 
#' @return Invisibly returns TRUE if successful, FALSE otherwise
#' 
#' @details
#' This function performs two steps:
#' \enumerate{
#'   \item Installs the cmdstanr R package from Stan's r-universe
#'   \item Installs the CmdStan C++ compiler (~200MB, 5-10 min compile time)
#' }
#' 
#' CmdStan compilation requires a C++ toolchain:
#' \itemize{
#'   \item macOS: Install Xcode Command Line Tools
#'   \item Windows: Install RTools
#'   \item Linux: gcc/g++ usually pre-installed
#' }
#' 
#' @examples
#' \dontrun{
#' # Install everything (recommended for first-time users)
#' install_bmm_dependencies()
#' 
#' # Non-interactive mode (for scripts)
#' install_bmm_dependencies(ask = FALSE)
#' 
#' # Only check, don't install
#' install_bmm_dependencies(install_cmdstan = FALSE)
#' }
#' 
#' @export
install_bmm_dependencies <- function(install_cmdstan = TRUE, 
                                     ask = TRUE,
                                     cores = parallel::detectCores(),
                                     verbose = TRUE) {
  
  if (verbose) {
    cat("\n╔════════════════════════════════════════════════════════╗\n")
    cat("║  BMM Dependency Installation                           ║\n")
    cat("╚════════════════════════════════════════════════════════╝\n\n")
  }
  
  all_ok <- TRUE
  
  # ============================================================================
  # Step 1: Check/Install cmdstanr R package
  # ============================================================================
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    if (verbose) {
      cat("Step 1: Installing cmdstanr R package...\n")
      cat("  Source: https://stan-dev.r-universe.dev\n\n")
    }
    
    if (ask) {
      response <- readline("Install cmdstanr? (y/n): ")
      if (tolower(response) != "y") {
        cat("Installation cancelled.\n")
        return(invisible(FALSE))
      }
    }
    
    tryCatch({
      install.packages(
        "cmdstanr",
        repos = c("https://stan-dev.r-universe.dev", getOption("repos")),
        type = "source"
      )
      if (verbose) cat("\n✓ cmdstanr installed successfully\n\n")
    }, error = function(e) {
      if (verbose) {
        cat("\n✗ Failed to install cmdstanr\n")
        cat("Error:", e$message, "\n\n")
        cat("Try manually:\n")
        cat("  install.packages('cmdstanr',\n")
        cat("    repos = c('https://stan-dev.r-universe.dev',\n")
        cat("              getOption('repos')))\n\n")
      }
      all_ok <<- FALSE
      return(invisible(FALSE))
    })
  } else {
    if (verbose) {
      cat("Step 1: cmdstanr R package\n")
      cat(sprintf("  ✓ Already installed (version %s)\n\n", 
                  packageVersion("cmdstanr")))
    }
  }
  
  # ============================================================================
  # Step 2: Check/Install CmdStan binary
  # ============================================================================
  if (install_cmdstan && all_ok) {
    
    # Check if already installed - CORRECT METHOD
    cmdstan_exists <- check_cmdstan_installed()
    
    if (cmdstan_exists) {
      if (verbose) {
        version <- tryCatch({
          cmdstanr::cmdstan_version()
        }, error = function(e) {
          "unknown"
        })
        
        path <- tryCatch({
          cmdstanr::cmdstan_path()
        }, error = function(e) {
          "unknown"
        })
        
        cat("Step 2: CmdStan binary\n")
        cat(sprintf("  ✓ Already installed (version %s)\n", version))
        cat(sprintf("  Location: %s\n\n", path))
      }
    } else {
      # Need to install
      if (verbose) {
        cat("Step 2: Installing CmdStan binary...\n")
        cat("  Size: ~200MB download\n")
        cat("  Time: 5-10 minutes (compilation)\n")
        cat(sprintf("  Cores: %d\n\n", cores))
      }
      
      if (ask) {
        # Get default path safely
        default_path <- tryCatch({
          cmdstanr::cmdstan_default_install_path()
        }, error = function(e) {
          "~/.cmdstan"
        })
        
        cat("CmdStan will be installed to:\n")
        cat(sprintf("  %s\n\n", default_path))
        response <- readline("Continue? (y/n): ")
        if (tolower(response) != "y") {
          cat("\nInstallation cancelled.\n")
          cat("You can install later with:\n")
          cat("  cmdstanr::install_cmdstan()\n\n")
          return(invisible(FALSE))
        }
      }
      
      if (verbose) cat("\nDownloading and compiling...\n")
      
      tryCatch({
        cmdstanr::install_cmdstan(cores = cores)
        
        if (verbose) {
          version <- tryCatch({
            cmdstanr::cmdstan_version()
          }, error = function(e) {
            "unknown"
          })
          
          path <- tryCatch({
            cmdstanr::cmdstan_path()
          }, error = function(e) {
            "unknown"
          })
          
          cat("\n✓ CmdStan installed successfully\n")
          cat(sprintf("  Version: %s\n", version))
          cat(sprintf("  Location: %s\n\n", path))
        }
      }, error = function(e) {
        if (verbose) {
          cat("\n✗ Failed to install CmdStan\n")
          cat("Error:", e$message, "\n\n")
          cat("Common issues:\n")
          cat("  • Missing C++ compiler (install Xcode/RTools/gcc)\n")
          cat("  • Insufficient disk space (need ~500MB)\n")
          cat("  • Network timeout (retry with: cmdstanr::install_cmdstan())\n\n")
        }
        all_ok <<- FALSE
        return(invisible(FALSE))
      })
    }
  }
  
  # ============================================================================
  # Final Summary
  # ============================================================================
  if (verbose) {
    cat("═══════════════════════════════════════════════════════\n")
    
    # Check final status
    cmdstan_ready <- check_cmdstan_installed()
    
    if (all_ok && (!install_cmdstan || cmdstan_ready)) {
      cat("✓ All dependencies ready!\n")
      cat("═══════════════════════════════════════════════════════\n\n")
      
      cat("You can now use BMM:\n\n")
      cat("  library(BMM)\n")
      cat("  Sys.setenv(TBB_CXX_TYPE = 'gcc')  # Already set\n\n")
      cat("  # Compile Stan models (once per session)\n")
      cat("  compiled <- compile_stan_models()\n\n")
      cat("  # Run analysis\n")
      cat("  results <- run_mr_model_comparison(\n")
      cat("    beta_X, beta_Y, se_X, se_Y,\n")
      cat("    compiled_models = compiled\n")
      cat("  )\n\n")
      
      return(invisible(TRUE))
      
    } else {
      cat("⚠️  Installation incomplete\n")
      cat("═══════════════════════════════════════════════════════\n\n")
      
      if (!requireNamespace("cmdstanr", quietly = TRUE)) {
        cat("Missing: cmdstanr R package\n")
        cat("  Install: install.packages('cmdstanr',\n")
        cat("            repos = c('https://stan-dev.r-universe.dev',\n")
        cat("                     getOption('repos')))\n\n")
      }
      
      if (install_cmdstan && !check_cmdstan_installed()) {
        cat("Missing: CmdStan binary\n")
        cat("  Install: cmdstanr::install_cmdstan()\n\n")
      }
      
      return(invisible(FALSE))
    }
  }
  
  invisible(all_ok)
}

#' Check BMM Environment
#' 
#' Check if all BMM requirements are met.
#' 
#' @param verbose Logical, print detailed output (default: TRUE)
#' 
#' @return Invisibly returns TRUE if all requirements met, FALSE otherwise
#' 
#' @export
setup_bmm_environment <- function(verbose = TRUE) {
  
  if (verbose) {
    cat("\n═══ BMM Environment Check ═══\n\n")
  }
  
  all_ok <- TRUE
  
  # TBB environment
  tbb <- Sys.getenv("TBB_CXX_TYPE")
  if (tbb == "") {
    if (verbose) cat("⚠️  TBB_CXX_TYPE not set, setting now...\n")
    Sys.setenv(TBB_CXX_TYPE = "gcc")
    tbb <- "gcc"
  }
  if (verbose) cat(sprintf("✓ TBB_CXX_TYPE = '%s'\n", tbb))
  
  # cmdstanr package
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    if (verbose) {
      cat("✗ cmdstanr package not installed\n")
      cat("  Run: install_bmm_dependencies()\n")
    }
    all_ok <- FALSE
  } else {
    if (verbose) {
      cat(sprintf("✓ cmdstanr package (version %s)\n", 
                  packageVersion("cmdstanr")))
    }
    
    # CmdStan binary - CORRECT CHECK
    if (!check_cmdstan_installed()) {
      if (verbose) {
        cat("✗ CmdStan binary not installed\n")
        cat("  Run: install_bmm_dependencies()\n")
      }
      all_ok <- FALSE
    } else {
      if (verbose) {
        version <- tryCatch({
          cmdstanr::cmdstan_version()
        }, error = function(e) {
          "unknown"
        })
        cat(sprintf("✓ CmdStan (version %s)\n", version))
      }
    }
  }
  
  # Plink (optional)
  plink_paths <- c("/usr/bin/plink", "/usr/local/bin/plink", 
                   "/opt/homebrew/bin/plink", "~/bin/plink")
  plink_found <- FALSE
  for (p in plink_paths) {
    p_expanded <- path.expand(p)
    if (file.exists(p_expanded)) {
      if (verbose) cat(sprintf("✓ Plink found: %s\n", p_expanded))
      plink_found <- TRUE
      break
    }
  }
  if (!plink_found && verbose) {
    cat("ℹ  Plink not found (optional, for harmonization)\n")
  }
  
  # Summary
  if (verbose) {
    cat("\n")
    if (all_ok) {
      cat("✓ Environment ready!\n\n")
    } else {
      cat("⚠️  Some requirements missing\n")
      cat("Run: install_bmm_dependencies()\n\n")
    }
  }
  
  invisible(all_ok)
}