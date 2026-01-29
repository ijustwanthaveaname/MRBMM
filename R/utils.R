################################################################################
# S3 Methods for BMM Results
################################################################################

#' Summary Method for BMM Results
#' 
#' @param object bmm_result object from run_mr_model_comparison
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the object
#' @export
summary.bmm_result <- function(object, ...) {
  cat("\n")
  cat("══════════════════════════════════════════════════════════\n")
  cat(" BMM Analysis Results\n")
  cat("══════════════════════════════════════════════════════════\n\n")
  
  # Model info
  cat(sprintf("Selected model: %s\n", object$selected_model))
  
  if (!is.null(object$final_result$priors$prior_selected_ivs)) {
    cat(sprintf("Number of SNPs used: %d\n", 
                length(object$final_result$priors$prior_selected_ivs)))
  }
  
  # Causal effect
  cat("\n─── Causal Effect (theta) ───\n")
  cat(sprintf("  Estimate: %.4f\n", object$final_result$theta$mean))
  cat(sprintf("  Std Error: %.4f\n", object$final_result$theta$se))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n", 
              object$final_result$theta$ci[1], 
              object$final_result$theta$ci[2]))
  cat(sprintf("  P-value: %.4e %s\n", 
              object$final_result$theta$p_value,
              ifelse(object$final_result$theta$p_value < 0.05, "***", "")))
  
  if (object$final_result$theta$p_value < 0.05) {
    direction <- ifelse(object$final_result$theta$mean > 0, "POSITIVE", "NEGATIVE")
    cat(sprintf("  → Significant %s causal effect\n", direction))
  } else {
    cat("  → No significant causal effect\n")
  }
  
  # Confounding bias (4-comp only)
  if (object$selected_model == "4-comp" && !is.null(object$final_result$delta)) {
    cat("\n─── Confounding Bias (delta = beta2 - theta) ───\n")
    cat(sprintf("  Estimate: %.4f\n", object$final_result$delta$mean))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
                object$final_result$delta$ci[1],
                object$final_result$delta$ci[2]))
    
    if (object$final_result$delta$ci[1] > 0) {
      cat("  → Significant POSITIVE confounding\n")
    } else if (object$final_result$delta$ci[2] < 0) {
      cat("  → Significant NEGATIVE confounding\n")
    } else {
      cat("  → No significant confounding\n")
    }
  }
  
  # Component weights
  cat("\n─── Component Weights (Posterior Mean) ───\n")
  n_comp <- object$final_result$n_components
  comp_names <- if (n_comp == 3) {
    c("Valid", "Valid+UHP", "NULL/UHP")
  } else {
    c("Valid", "Valid+UHP", "Valid+UHP+CHP", "NULL")
  }
  
  for (i in 1:n_comp) {
    cat(sprintf("  π[%d] (%15s): %.3f (SD = %.3f)\n", 
                i, comp_names[i], 
                object$final_result$pi_posterior$mean[i],
                object$final_result$pi_posterior$sd[i]))
  }
  
  # Convergence
  cat("\n─── Convergence Diagnostics ───\n")
  if (object$final_result$convergence$converged) {
    cat("  ✓ All chains converged successfully\n")
    cat(sprintf("    Min ESS ratio: %.3f\n", 
                object$final_result$convergence$min_ess_ratio))
  } else {
    cat("  ⚠️  Convergence issues detected:\n")
    for (issue in object$final_result$convergence$issues) {
      cat(sprintf("     - %s\n", issue))
    }
  }
  
  # Timing
  if (!is.null(object$elapsed_time_minutes)) {
    cat("\n─── Timing ───\n")
    cat(sprintf("  Total time: %.2f minutes\n", object$elapsed_time_minutes))
    if (!is.null(object$parallel_used) && object$parallel_used) {
      cat("  Parallel: Yes\n")
    }
  }
  
  cat("\n══════════════════════════════════════════════════════════\n\n")
  
  invisible(object)
}

#' Print Method for BMM Results
#' 
#' @param x bmm_result object
#' @param ... Additional arguments passed to summary
#' @return Invisibly returns x
#' @export
print.bmm_result <- function(x, ...) {
  summary(x, ...)
}

#' Extract theta Estimates
#' 
#' @param object bmm_result object
#' @return Data frame with theta estimates from both models
#' @export
get_theta_estimates <- function(object) {
  data.frame(
    model = c("3-comp", "4-comp", "Selected"),
    theta = c(
      object$result_3comp$theta$mean,
      object$result_4comp$theta$mean,
      object$final_result$theta$mean
    ),
    se = c(
      object$result_3comp$theta$se,
      object$result_4comp$theta$se,
      object$final_result$theta$se
    ),
    ci_lower = c(
      object$result_3comp$theta$ci[1],
      object$result_4comp$theta$ci[1],
      object$final_result$theta$ci[1]
    ),
    ci_upper = c(
      object$result_3comp$theta$ci[2],
      object$result_4comp$theta$ci[2],
      object$final_result$theta$ci[2]
    ),
    p_value = c(
      object$result_3comp$theta$p_value,
      object$result_4comp$theta$p_value,
      object$final_result$theta$p_value
    )
  )
}
