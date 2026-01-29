#' Weighted mode estimate
#'
#' Internal weighted mode estimator used in BMM.
#'
#' @keywords internal
weighted_mode_estimate <- function(beta_X, beta_Y, se_X, se_Y, 
                                   phi = 1, 
                                   assume_NOME = FALSE,
                                   verbose = TRUE) {
  
  if (verbose) cat("\n=== Weighted Mode Estimation (Hartwig 2017) ===\n")
  
  # Data cleaning
  valid_idx <- is.finite(beta_X) & is.finite(beta_Y) & 
               se_X > 0 & se_Y > 0 & abs(beta_X) > 1e-6
  
  if (sum(valid_idx) < length(beta_X)) {
    n_removed <- length(beta_X) - sum(valid_idx)
    if (verbose) cat(sprintf("  Removed %d SNPs with invalid data\n", n_removed))
  }
  
  beta_X <- beta_X[valid_idx]
  beta_Y <- beta_Y[valid_idx]
  se_X <- se_X[valid_idx]
  se_Y <- se_Y[valid_idx]
  J <- length(beta_X)
  
  if (J < 3) {
    if (verbose) cat("  ⚠️  Insufficient SNPs (< 3)\n")
    return(list(
      theta_mode = NA, 
      se_mode = NA,
      ci_lower = NA,
      ci_upper = NA,
      converged = FALSE,
      n_snps = J
    ))
  }
  
  # === Ratio estimates (Eq. 3 in Hartwig 2017) ===
  ratio_estimates <- beta_Y / beta_X
  
  # === Standard errors (Eq. 4 in Hartwig 2017) ===
  if (assume_NOME) {
    # Simplified under NOME assumption: σ_Rj = σ_Yj / |β̂_Xj|
    se_ratio <- se_Y / abs(beta_X)
  } else {
    # Full delta method formula (recommended)
    se_ratio <- sqrt(
      se_Y^2 / beta_X^2 + 
      beta_Y^2 * se_X^2 / beta_X^4
    )
  }
  
  # === Weights (Eq. 5 - inverse variance weighting) ===
  weights <- se_ratio^(-2)
  weights <- weights / sum(weights)  # Normalize to sum to 1
  
  # === Search grid for mode ===
  q05 <- stats::quantile(ratio_estimates, 0.05, na.rm = TRUE)
  q95 <- stats::quantile(ratio_estimates, 0.95, na.rm = TRUE)
  extend <- (q95 - q05) * 0.5
  theta_grid <- seq(q05 - extend, q95 + extend, length.out = 100)
  
  # === Bandwidth selection (Eq. 7 - Modified Silverman's rule) ===
  sd_ratios <- stats::sd(ratio_estimates, na.rm = TRUE)
  mad_ratios <- stats::median(abs(ratio_estimates - stats::median(ratio_estimates, na.rm = TRUE)), 
                       na.rm = TRUE)
  
  # s = 0.9 * min(sd, 1.4826*MAD) / L^(1/5)
  s <- 0.9 * min(sd_ratios, 1.4826 * mad_ratios) * J^(-1/5)
  bandwidth <- phi * s
  
  if (verbose) {
    cat(sprintf("  Number of SNPs: %d\n", J))
    cat(sprintf("  Ratio estimates: median = %.4f, IQR = [%.4f, %.4f]\n",
                stats::median(ratio_estimates, na.rm = TRUE),
                stats::quantile(ratio_estimates, 0.25, na.rm = TRUE),
                stats::quantile(ratio_estimates, 0.75, na.rm = TRUE)))
    cat(sprintf("  Bandwidth: h = %.4f (φ = %.2f, s = %.4f)\n", bandwidth, phi, s))
    cat(sprintf("  Search range: [%.4f, %.4f]\n", min(theta_grid), max(theta_grid)))
  }
  
  # === Kernel density estimation (Eq. 6 in Hartwig 2017) ===
  # f(x) = (1 / h√(2π)) Σ wj·exp[-1/2 ((x - β̂_Rj) / h)²]
  density_values <- sapply(theta_grid, function(theta) {
    kernel_vals <- stats::dnorm((ratio_estimates - theta) / bandwidth)
    sum(weights * kernel_vals)
  })
  
  # === Find mode (maximum density point) ===
  max_idx <- which.max(density_values)
  theta_mode <- theta_grid[max_idx]
  max_density <- density_values[max_idx]
  
  # === Confidence interval (based on density threshold method) ===
  # Use 50% of max density as threshold
  threshold_density <- max_density * 0.5
  in_ci <- density_values > threshold_density
  
  if (sum(in_ci) > 0) {
    theta_ci_lower <- min(theta_grid[in_ci])
    theta_ci_upper <- max(theta_grid[in_ci])
  } else {
    # Fallback: use ±1.96*bandwidth
    theta_ci_lower <- theta_mode - 1.96 * bandwidth
    theta_ci_upper <- theta_mode + 1.96 * bandwidth
  }
  
  # Approximate standard error from CI width
  se_mode <- (theta_ci_upper - theta_ci_lower) / (2 * 1.96)
  
  if (verbose) {
    cat(sprintf("  ✓ Mode estimate: %.4f\n", theta_mode))
    cat(sprintf("  ✓ Approx SE: %.4f\n", se_mode))
    cat(sprintf("  ✓ 95%% CI: [%.4f, %.4f]\n", theta_ci_lower, theta_ci_upper))
  }
  
  return(list(
    theta_mode = theta_mode,
    se_mode = se_mode,
    ci_lower = theta_ci_lower,
    ci_upper = theta_ci_upper,
    bandwidth = bandwidth,
    phi = phi,
    theta_grid = theta_grid,
    density_values = density_values,
    ratio_estimates = ratio_estimates,
    weights = weights,
    n_snps = J,
    assume_NOME = assume_NOME,
    converged = TRUE,
    method = "Weighted Mode (Hartwig 2017)"
  ))
}

################################################################################
# PART 2: model definition
################################################################################

# 3comp
stan_model_3comp <- "
data {
  int<lower=0> J;
  vector[J] beta_X_hat;
  vector[J] beta_Y_hat;
  vector<lower=0>[J] se_X;
  vector<lower=0>[J] se_Y;
  real<lower=-1,upper=1> rho_e;
  
  real<lower=0> sigma_theta_c;
  real<lower=0> a_x;
  real<lower=0> b_x;
  real<lower=0> a_y;
  real<lower=0> b_y;
  vector<lower=0>[3] alpha_pi;
}

parameters {
  real theta_raw;
  real<lower=0> sigma_x_sq;
  real<lower=0> sigma_y_sq;
  simplex[3] pi;
}

transformed parameters {
  real theta = theta_raw * sigma_theta_c;
}

model {
  real epsilon = 1e-6;
  array[3] matrix[2,2] V;
  matrix[2,2] Sigma_obs;
  vector[2] beta_hat;
  vector[3] log_lik_components;
  
  theta_raw ~ std_normal();
  sigma_x_sq ~ inv_gamma(a_x, b_x);
  sigma_y_sq ~ inv_gamma(a_y, b_y);
  pi ~ dirichlet(alpha_pi);
  
  V[1][1,1] = sigma_x_sq;
  V[1][1,2] = theta * sigma_x_sq;
  V[1][2,1] = theta * sigma_x_sq;
  V[1][2,2] = theta * theta * sigma_x_sq;
  
  V[2][1,1] = sigma_x_sq;
  V[2][1,2] = theta * sigma_x_sq;
  V[2][2,1] = theta * sigma_x_sq;
  V[2][2,2] = theta * theta * sigma_x_sq + sigma_y_sq;
  
  V[3][1,1] = epsilon;
  V[3][1,2] = 0;
  V[3][2,1] = 0;
  V[3][2,2] = sigma_y_sq;
  
  for (j in 1:J) {
    beta_hat[1] = beta_X_hat[j];
    beta_hat[2] = beta_Y_hat[j];
    
    Sigma_obs[1,1] = se_X[j]^2;
    Sigma_obs[1,2] = rho_e * se_X[j] * se_Y[j];
    Sigma_obs[2,1] = rho_e * se_X[j] * se_Y[j];
    Sigma_obs[2,2] = se_Y[j]^2;
    
    for (k in 1:3) {
      log_lik_components[k] = multi_normal_lpdf(beta_hat | [0,0]', V[k] + Sigma_obs);
    }
    
    target += log_sum_exp(log(pi) + log_lik_components);
  }
}

generated quantities {
  vector[J] log_lik;
  
  for (j in 1:J) {
    vector[2] beta_hat;
    matrix[2,2] Sigma_obs;
    vector[3] log_lik_components;
    array[3] matrix[2,2] V;
    real epsilon = 1e-6;
    
    beta_hat[1] = beta_X_hat[j];
    beta_hat[2] = beta_Y_hat[j];
    
    Sigma_obs[1,1] = se_X[j]^2;
    Sigma_obs[1,2] = rho_e * se_X[j] * se_Y[j];
    Sigma_obs[2,1] = rho_e * se_X[j] * se_Y[j];
    Sigma_obs[2,2] = se_Y[j]^2;
    
    V[1][1,1] = sigma_x_sq;
    V[1][1,2] = theta * sigma_x_sq;
    V[1][2,1] = theta * sigma_x_sq;
    V[1][2,2] = theta * theta * sigma_x_sq;
    
    V[2][1,1] = sigma_x_sq;
    V[2][1,2] = theta * sigma_x_sq;
    V[2][2,1] = theta * sigma_x_sq;
    V[2][2,2] = theta * theta * sigma_x_sq + sigma_y_sq;
    
    V[3][1,1] = epsilon;
    V[3][1,2] = 0;
    V[3][2,1] = 0;
    V[3][2,2] = sigma_y_sq;
    
    for (k in 1:3) {
      log_lik_components[k] = multi_normal_lpdf(beta_hat | [0,0]', V[k] + Sigma_obs);
    }
    
    log_lik[j] = log_sum_exp(log(pi) + log_lik_components);
  }
}
"

# 4comp
stan_model_4comp <- "
data {
  int<lower=0> J;
  vector[J] beta_X_hat;
  vector[J] beta_Y_hat;
  vector<lower=0>[J] se_X;
  vector<lower=0>[J] se_Y;
  real<lower=-1,upper=1> rho_e;
  
  real<lower=0> sigma_theta_c;
  real<lower=0> a_x;
  real<lower=0> b_x;
  real<lower=0> a_y;
  real<lower=0> b_y;
  real<lower=0> a_tilde_c;
  real<lower=0> b_tilde_c;
  vector<lower=0>[4] alpha_pi;
}

parameters {
  real theta_raw;
  real beta2_raw;
  real<lower=0> sigma_x_sq;
  real<lower=0> sigma_y_sq;
  real<lower=0> sigma_tilde_c_sq;
  simplex[4] pi;
}

transformed parameters {
  real theta = theta_raw * sigma_theta_c;
  real beta2 = beta2_raw * sigma_theta_c;
  real delta = beta2 - theta;
  
  real<lower=0> tau1_sq = sigma_y_sq;
  real<lower=tau1_sq> tau2_sq = sigma_y_sq + sigma_tilde_c_sq;
}

model {
  real epsilon = 1e-6;
  array[4] matrix[2,2] V;
  matrix[2,2] Sigma_obs;
  vector[2] beta_hat;
  vector[4] log_lik_components;
  
  theta_raw ~ std_normal();
  beta2_raw ~ std_normal();
  sigma_x_sq ~ inv_gamma(a_x, b_x);
  sigma_y_sq ~ inv_gamma(a_y, b_y);
  sigma_tilde_c_sq ~ inv_gamma(a_tilde_c, b_tilde_c);
  pi ~ dirichlet(alpha_pi);
  
  V[1][1,1] = sigma_x_sq;
  V[1][1,2] = theta * sigma_x_sq;
  V[1][2,1] = theta * sigma_x_sq;
  V[1][2,2] = theta * theta * sigma_x_sq;
  
  V[2][1,1] = sigma_x_sq;
  V[2][1,2] = theta * sigma_x_sq;
  V[2][2,1] = theta * sigma_x_sq;
  V[2][2,2] = theta * theta * sigma_x_sq + sigma_y_sq;
  
  V[3][1,1] = sigma_x_sq;
  V[3][1,2] = beta2 * sigma_x_sq;
  V[3][2,1] = beta2 * sigma_x_sq;
  V[3][2,2] = beta2 * beta2 * sigma_x_sq + tau2_sq;
  
  V[4][1,1] = epsilon;
  V[4][1,2] = 0;
  V[4][2,1] = 0;
  V[4][2,2] = sigma_y_sq;
  
  for (j in 1:J) {
    beta_hat[1] = beta_X_hat[j];
    beta_hat[2] = beta_Y_hat[j];
    
    Sigma_obs[1,1] = se_X[j]^2;
    Sigma_obs[1,2] = rho_e * se_X[j] * se_Y[j];
    Sigma_obs[2,1] = rho_e * se_X[j] * se_Y[j];
    Sigma_obs[2,2] = se_Y[j]^2;
    
    for (k in 1:4) {
      log_lik_components[k] = multi_normal_lpdf(beta_hat | [0,0]', V[k] + Sigma_obs);
    }
    
    target += log_sum_exp(log(pi) + log_lik_components);
  }
}

generated quantities {
  vector[J] log_lik;
  
  for (j in 1:J) {
    vector[2] beta_hat;
    matrix[2,2] Sigma_obs;
    vector[4] log_lik_components;
    array[4] matrix[2,2] V;
    real epsilon = 1e-6;
    
    beta_hat[1] = beta_X_hat[j];
    beta_hat[2] = beta_Y_hat[j];
    
    Sigma_obs[1,1] = se_X[j]^2;
    Sigma_obs[1,2] = rho_e * se_X[j] * se_Y[j];
    Sigma_obs[2,1] = rho_e * se_X[j] * se_Y[j];
    Sigma_obs[2,2] = se_Y[j]^2;
    
    V[1][1,1] = sigma_x_sq;
    V[1][1,2] = theta * sigma_x_sq;
    V[1][2,1] = theta * sigma_x_sq;
    V[1][2,2] = theta * theta * sigma_x_sq;
    
    V[2][1,1] = sigma_x_sq;
    V[2][1,2] = theta * sigma_x_sq;
    V[2][2,1] = theta * sigma_x_sq;
    V[2][2,2] = theta * theta * sigma_x_sq + sigma_y_sq;
    
    V[3][1,1] = sigma_x_sq;
    V[3][1,2] = beta2 * sigma_x_sq;
    V[3][2,1] = beta2 * sigma_x_sq;
    V[3][2,2] = beta2 * beta2 * sigma_x_sq + sigma_y_sq + sigma_tilde_c_sq;
    
    V[4][1,1] = epsilon;
    V[4][1,2] = 0;
    V[4][2,1] = 0;
    V[4][2,2] = sigma_y_sq;
    
    for (k in 1:4) {
      log_lik_components[k] = multi_normal_lpdf(beta_hat | [0,0]', V[k] + Sigma_obs);
    }
    
    log_lik[j] = log_sum_exp(log(pi) + log_lik_components);
  }
}
"

################################################################################
# PART 3: compile stan models
################################################################################

compile_stan_models <- function(force_recompile = FALSE) {
  cat("\n===== Compiling Stan Models =====\n")
  
  stan_dir <- tempdir()
  
  model_3comp_file <- file.path(stan_dir, "mr_3comp.stan")
  writeLines(stan_model_3comp, model_3comp_file)
  cat("Compiling 3-component model...\n")
  model_3comp <- cmdstanr::cmdstan_model(model_3comp_file, force_recompile = force_recompile)
  
  model_4comp_file <- file.path(stan_dir, "mr_4comp.stan")
  writeLines(stan_model_4comp, model_4comp_file)
  cat("Compiling 4-component model...\n")
  model_4comp <- cmdstanr::cmdstan_model(model_4comp_file, force_recompile = force_recompile)
  
  cat("Compilation complete!\n")
  
  return(list(
    model_3comp = model_3comp,
    model_4comp = model_4comp
  ))
}

################################################################################
# PART 4: LOO calculation
################################################################################

compute_robust_loo <- function(fit, log_lik_matrix, k_threshold = 0.7, 
                               reloo_threshold = 0.1) {
  
  cat("\nComputing LOO-CV...\n")
  
  loo_result <- loo::loo(log_lik_matrix)
  k_vals <- loo_result$diagnostics$pareto_k
  n_high_k <- sum(k_vals > k_threshold)
  prop_high_k <- n_high_k / length(k_vals)
  
  if (n_high_k > 0) {
    cat(sprintf("  Warning: %d observations (%.1f%%) with Pareto k > %.1f\n", 
                n_high_k, 100*prop_high_k, k_threshold))
    
    if (prop_high_k > reloo_threshold) {
      cat("\n  → Applying moment matching...\n")
      
      tryCatch({
        loo_result <- loo::loo_moment_match(fit, loo = loo_result, 
                                       unconstrain_pars = TRUE)
        
        k_vals_new <- loo_result$diagnostics$pareto_k
        n_high_k_new <- sum(k_vals_new > k_threshold)
        
        if (n_high_k_new < n_high_k) {
          cat(sprintf("    ✓ Improved: %d → %d problematic observations\n", 
                      n_high_k, n_high_k_new))
        }
        
      }, error = function(e) {
        cat("    Warning: Moment matching failed\n")
        cat("    → Use results with caution\n")
      })
      
    } else {
      cat("    Note: Impact is minor (<10%), continuing\n")
    }
  } else {
    cat("    ✓ All Pareto k values < 0.7 (excellent)\n")
  }
  
  cat(sprintf("  ELPD = %.2f (SE = %.2f)\n", 
              loo_result$estimates["elpd_loo", "Estimate"],
              loo_result$estimates["elpd_loo", "SE"]))
  
  return(loo_result)
}

################################################################################
# PART 5: Adaptive thresholds
################################################################################
determine_adaptive_thresholds <- function(p_X, p_Y, 
                                          user_pval_threshold_X = NULL,
                                          user_pval_threshold_Y = NULL) {
    
    cat("\n===== Adaptive Threshold Determination =====\n")
    J <- length(p_X)
    
    if (!is.null(user_pval_threshold_X)) {
        thresh_X_strong <- user_pval_threshold_X
        cat(sprintf("Using user X threshold: %.2e\n", thresh_X_strong))
    } else {
        gwas_thresholds <- c(5e-8, 1e-6, 1e-5, 1e-4, 1e-3, 0.01, 0.05)
        min_ivs <- 10
        
        n_below <- sapply(gwas_thresholds, function(t) sum(p_X < t))
        valid_thresholds_idx <- which(n_below >= min_ivs)
        
        if (length(valid_thresholds_idx) > 0) {
            thresh_X_strong <- gwas_thresholds[valid_thresholds_idx][1] 
            cat(sprintf("Auto X threshold: %.2e (n=%d SNPs)\n", 
                        thresh_X_strong, sum(p_X < thresh_X_strong)))
        } else {
            thresh_X_strong <- stats::quantile(p_X, 0.10)
            cat(sprintf("Using 10th percentile: %.2e (n=%d SNPs)\n", 
                        thresh_X_strong, sum(p_X < thresh_X_strong)))
        }
    }
    
    if (!is.null(user_pval_threshold_Y)) {
        thresh_Y_sig <- user_pval_threshold_Y
        cat(sprintf("Using user Y threshold: %.2e\n", thresh_Y_sig))
    } else {
        if (J > 0) {
            thresh_Y_sig_bonf <- 0.05 / J
            thresh_Y_sig <- min(0.05, thresh_Y_sig_bonf)
            cat(sprintf("Auto Y threshold: %.2e (Bonferroni/0.05)\n", thresh_Y_sig))
        } else {
            thresh_Y_sig <- 0.05
            cat(sprintf("Using default Y threshold: %.2e\n", thresh_Y_sig))
        }
    }
    
    thresh_weak <- 0.05 
    cat(sprintf("Weak X threshold: %.2e\n", thresh_weak))
    
    return(list(
        thresh_X_strong = thresh_X_strong,
        thresh_Y_sig = thresh_Y_sig,
        thresh_weak = thresh_weak
    ))
}

################################################################################
# PART 6: Dirichlet prior
################################################################################

compute_smart_dirichlet_prior <- function(beta_X, beta_Y, se_X, se_Y, 
                                         p_X = NULL, p_Y = NULL,
                                         pval_threshold_X = NULL,
                                         pval_threshold_Y = NULL,
                                         prior_strength = 1.0,
                                         n_components = 3) {
  
  J <- length(beta_X)
  
  if (is.null(p_X)) p_X <- 2 * stats::pnorm(-abs(beta_X / se_X))
  if (is.null(p_Y)) p_Y <- 2 * stats::pnorm(-abs(beta_Y / se_Y))
  
  cat(sprintf("\n===== Smart Dirichlet Prior (%d components) =====\n", n_components))
  
  thresholds <- determine_adaptive_thresholds(p_X, p_Y, pval_threshold_X, pval_threshold_Y)
  thresh_X_strong <- thresholds$thresh_X_strong
  thresh_Y_sig <- thresholds$thresh_Y_sig
  thresh_weak <- thresholds$thresh_weak
  
  X_strong <- p_X < thresh_X_strong
  X_weak <- p_X > thresh_weak
  Y_sig <- p_Y < thresh_Y_sig
  Y_weak <- p_Y > thresh_weak
  
  n_quad1 <- sum(X_strong & Y_sig)
  n_quad2 <- sum(X_strong & Y_weak)
  n_quad3 <- sum(X_weak & Y_sig)
  n_quad4 <- sum(X_weak & Y_weak)
  
  cat(sprintf("\nQuadrant classification (J=%d):\n", J))
  cat(sprintf("  ① X-strong & Y-sig:  %3d (%.1f%%)\n", n_quad1, 100*n_quad1/J))
  cat(sprintf("  ② X-strong & Y-weak: %3d (%.1f%%)\n", n_quad2, 100*n_quad2/J))
  cat(sprintf("  ③ X-weak & Y-sig:    %3d (%.1f%%)\n", n_quad3, 100*n_quad3/J))
  cat(sprintf("  ④ X-weak & Y-weak:   %3d (%.1f%%)\n", n_quad4, 100*n_quad4/J))
  
  X_strong_idx <- which(X_strong)
  pleiotropy_detected <- FALSE
  ratio_est <- 0
  
  if (length(X_strong_idx) >= 5) {
    valid_idx <- X_strong_idx[is.finite(beta_X[X_strong_idx]) & 
                              is.finite(beta_Y[X_strong_idx]) &
                              abs(beta_X[X_strong_idx]) > 1e-10]
    
    if (length(valid_idx) >= 5) {
      weights <- 1 / (se_Y[valid_idx]^2)
      ratio_est <- sum(weights * beta_Y[valid_idx] * beta_X[valid_idx]) / 
                   sum(weights * beta_X[valid_idx]^2)
      
      X_mat <- cbind(1, beta_X[valid_idx])
      W <- diag(weights)
      tryCatch({
        coef_egger <- solve(t(X_mat) %*% W %*% X_mat) %*% t(X_mat) %*% W %*% beta_Y[valid_idx]
        residuals <- beta_Y[valid_idx] - X_mat %*% coef_egger
        mse <- sum(weights * residuals^2) / (length(valid_idx) - 2)
        var_coef <- mse * solve(t(X_mat) %*% W %*% X_mat)
        se_intercept <- sqrt(var_coef[1, 1])
        z_intercept <- abs(coef_egger[1] / se_intercept)
        
        pleiotropy_detected <- (z_intercept > 1.96)
        
        cat(sprintf("\nMR estimates (n=%d IVs):\n", length(valid_idx)))
        cat(sprintf("  IVW ratio: %.4f\n", ratio_est))
        cat(sprintf("  Egger intercept: %.4f (Z=%.2f)\n", coef_egger[1], z_intercept))
        cat(sprintf("  Pleiotropy: %s\n", ifelse(pleiotropy_detected, "Detected", "Not detected")))
      }, error = function(e) {
        cat(sprintf("\nIVW ratio: %.4f\n", ratio_est))
      })
    }
  }
  
  if (n_components == 3) {
    alpha_pi <- rep(1, 3)
    comp_names <- c("Valid", "Valid+UHP", "NULL/UHP")
    
    cat(sprintf("\nDirichlet heuristics (3-comp):\n"))
    
    prop_null <- n_quad4 / J
    if (prop_null > 0.2) {
      alpha_pi[3] <- alpha_pi[3] + 3 * (prop_null / 0.2)
      cat(sprintf("  [1] Null %.1f%% → α[NULL/UHP] += %.2f\n",
                  100*prop_null, 3*(prop_null/0.2)))
    }
    
    prop_uhp <- n_quad3 / J
    if (prop_uhp > 0.05) {
      alpha_pi[3] <- alpha_pi[3] + 2 * (prop_uhp / 0.05)
      cat(sprintf("  [2] Weak-X strong-Y %.1f%% → α[NULL/UHP] += %.2f\n",
                  100*prop_uhp, 2*(prop_uhp/0.05)))
    }
    
    prop_ambiguous <- n_quad1 / J
    if (prop_ambiguous > 0.1) {
      if (pleiotropy_detected) {
        alpha_pi[2] <- alpha_pi[2] + 3 * (prop_ambiguous / 0.1)
        cat(sprintf("  [3] X-strong Y-strong %.1f%% with pleiotropy → α[Valid+UHP] += %.2f\n",
                    100*prop_ambiguous, 3*(prop_ambiguous/0.1)))
      } else {
        alpha_pi[1] <- alpha_pi[1] + 2 * (prop_ambiguous / 0.1)
        alpha_pi[2] <- alpha_pi[2] + 1 * (prop_ambiguous / 0.1)
        cat(sprintf("  [3] X-strong Y-strong %.1f%% no pleiotropy → α[Valid] += 2, α[Valid+UHP] += 1\n",
                    100*prop_ambiguous))
      }
    }
    
    prop_valid_clean <- n_quad2 / J
    if (prop_valid_clean > 0.1) {
      alpha_pi[1] <- alpha_pi[1] + 4 * (prop_valid_clean / 0.1)
      cat(sprintf("  [4] X-strong Y-weak %.1f%% → α[Valid] += %.2f\n",
                  100*prop_valid_clean, 4*(prop_valid_clean/0.1)))
    }
    
  } else if (n_components == 4) {
    alpha_pi <- rep(1, 4)
    comp_names <- c("Valid", "Valid+UHP", "Valid+UHP+CHP", "NULL")
    
    cat(sprintf("\nDirichlet heuristics (4-comp):\n"))
    
    prop_null <- n_quad4 / J
    if (prop_null > 0.2) {
      alpha_pi[4] <- alpha_pi[4] + 3 * (prop_null / 0.2)
      cat(sprintf("  [1] Null %.1f%% → α[NULL] += %.2f\n",
                  100*prop_null, 3*(prop_null/0.2)))
    }
    
    prop_uhp <- n_quad3 / J
    if (prop_uhp > 0.05) {
      if (pleiotropy_detected) {
        alpha_pi[2] <- alpha_pi[2] + 2 * (prop_uhp / 0.05)
        cat(sprintf("  [2] Weak-X strong-Y %.1f%% with pleiotropy → α[Valid+UHP] += %.2f\n",
                    100*prop_uhp, 2*(prop_uhp/0.05)))
      } else {
        alpha_pi[3] <- alpha_pi[3] + 2 * (prop_uhp / 0.05)
        cat(sprintf("  [2] Weak-X strong-Y %.1f%% no pleiotropy → α[UHP+CHP] += %.2f\n",
                    100*prop_uhp, 2*(prop_uhp/0.05)))
      }
    }
    
    prop_valid_clean <- n_quad2 / J
    if (prop_valid_clean > 0.1) {
      alpha_pi[1] <- alpha_pi[1] + 5 * (prop_valid_clean / 0.1)
      cat(sprintf("  [3] Clean Valid IVs %.1f%% → α[Valid] += %.2f\n",
                  100*prop_valid_clean, 5*(prop_valid_clean/0.1)))
    }
    
    prop_ambiguous <- n_quad1 / J
    if (prop_ambiguous > 0.1) {
      if (pleiotropy_detected) {
        alpha_pi[2] <- alpha_pi[2] + 3 * (prop_ambiguous / 0.1)
        cat(sprintf("  [4] X-strong Y-strong %.1f%% with pleiotropy → α[Valid+UHP] += %.2f\n",
                    100*prop_ambiguous, 3*(prop_ambiguous/0.1)))
      } else {
        alpha_pi[1] <- alpha_pi[1] + 3 * (prop_ambiguous / 0.1)
        alpha_pi[3] <- alpha_pi[3] + 3 * (prop_ambiguous / 0.1)
        cat(sprintf("  [4] X-strong Y-strong %.1f%% no pleiotropy → α[Valid] and α[UHP+CHP] += 3\n",
                    100*prop_ambiguous))
      }
    }
    
    p_ratio <- stats::median(p_Y) / stats::median(p_X)
    if (p_ratio > 10 && n_quad2 > n_quad1) {
      alpha_pi[1] <- alpha_pi[1] + 2
      cat(sprintf("  [5] Y much weaker (ratio=%.1f) → α[Valid] += 2\n", p_ratio))
    }
    
    if (p_ratio < 0.5 && n_quad3 > 0.1 * J) {
      alpha_pi[3] <- alpha_pi[3] + 2
      cat(sprintf("  [6] Y stronger (ratio=%.2f) + weak-X strong-Y → α[UHP+CHP] += 2\n",
                  p_ratio))
    }
  }
  
  if (J < 50) {
    alpha_pi <- alpha_pi * 2
    cat(sprintf("  [*] Small sample (J=%d) → all α × 2\n", J))
  } else if (J > 500) {
    alpha_pi <- alpha_pi * 0.5
    cat(sprintf("  [*] Large sample (J=%d) → all α × 0.5\n", J))
  }
  
  alpha_pi <- alpha_pi * prior_strength
  if (prior_strength != 1.0) {
    cat(sprintf("  [*] User strength: × %.2f\n", prior_strength))
  }
  
  alpha_pi <- pmax(alpha_pi, 0.5)
  
  cat(sprintf("\nFinal alpha_pi:\n"))
  for (i in 1:n_components) {
    cat(sprintf("  Component %d [%15s]: α=%.2f → E[π]=%.3f\n",
                i, comp_names[i], alpha_pi[i], alpha_pi[i]/sum(alpha_pi)))
  }
  
  prior_concentration <- sum(alpha_pi)
  cat(sprintf("\nPrior concentration: %.2f ", prior_concentration))
  if (prior_concentration < 8) {
    cat("(weak)\n")
  } else if (prior_concentration < 20) {
    cat("(moderate)\n")
  } else {
    cat("(strong)\n")
  }
  
  return(list(
    alpha_pi = alpha_pi,
    thresholds = thresholds,
    classification = list(
      n_quad1 = n_quad1, n_quad2 = n_quad2,
      n_quad3 = n_quad3, n_quad4 = n_quad4
    ),
    p_ratio = if (exists("p_ratio")) p_ratio else stats::median(p_Y) / stats::median(p_X),
    ratio_est = ratio_est,
    pleiotropy_detected = pleiotropy_detected,
    prior_concentration = prior_concentration
  ))
}


################################################################################
# PART 7: Empirical priors
################################################################################

compute_empirical_priors <- function(beta_X, beta_Y, se_X, se_Y, 
                                     p_X = NULL, 
                                     n_components = 3,
                                     use_weighted_mode = TRUE,
                                     min_ivs_for_mode = 10) {
    
    J <- length(beta_X)
    if (is.null(p_X)) p_X <- 2 * stats::pnorm(-abs(beta_X / se_X))
    
    cat("\n===== Computing Empirical Priors (Dual IV Selection Strategy) =====\n")
    
    # ========== 策略1: 严格筛选IV用于Weighted Mode ==========
    theta_robust <- NA
    mode_result <- NULL
    mode_selected_ivs <- NULL
    mode_threshold <- NA
    
    if (use_weighted_mode) {
        cat("\n[Strategy 1] Strict IV Selection for Weighted Mode (theta initialization):\n")
        
        # 严格的阶梯式p值阈值
        strict_thresholds <- c(5e-8, 1e-8, 5e-7, 1e-7, 5e-6, 1e-6, 5e-5, 1e-5, 5e-4, 1e-4, 5e-3, 1e-3)
        
        for (threshold in strict_thresholds) {
            potential_ivs <- which(p_X < threshold)
            n_potential <- length(potential_ivs)
            
            cat(sprintf("  Trying p_X < %.2e: %d IVs\n", threshold, n_potential))
            
            if (n_potential >= min_ivs_for_mode) {
                mode_selected_ivs <- potential_ivs
                mode_threshold <- threshold
                cat(sprintf("  ✓ Selected for Mode: %.2e (%d IVs)\n", 
                           mode_threshold, length(mode_selected_ivs)))
                break
            }
        }
        
        # 如果严格阈值都不满足，使用top N个SNP
        if (is.null(mode_selected_ivs)) {
            n_to_select <- min(min_ivs_for_mode, J)
            mode_selected_ivs <- order(p_X)[1:n_to_select]
            mode_threshold <- max(p_X[mode_selected_ivs])
            cat(sprintf("  ⚠️  Using top %d SNPs (p_X < %.2e)\n", 
                       n_to_select, mode_threshold))
        }
        
        # 质量控制
        valid_mode_ivs <- mode_selected_ivs[
            abs(beta_X[mode_selected_ivs]) > 1e-6 & 
            is.finite(beta_X[mode_selected_ivs]) & 
            is.finite(beta_Y[mode_selected_ivs])
        ]
        
        if (length(valid_mode_ivs) < length(mode_selected_ivs)) {
            cat(sprintf("  → QC: Removed %d SNPs with |beta_X| ≈ 0 or NA\n", 
                       length(mode_selected_ivs) - length(valid_mode_ivs)))
        }
        
        mode_selected_ivs <- valid_mode_ivs
        
        # 计算Weighted Mode
        if (length(mode_selected_ivs) >= 3) {
            mode_result <- weighted_mode_estimate(
                beta_X = beta_X[mode_selected_ivs], 
                beta_Y = beta_Y[mode_selected_ivs], 
                se_X = se_X[mode_selected_ivs], 
                se_Y = se_Y[mode_selected_ivs], 
                phi = 1, 
                assume_NOME = FALSE,
                verbose = FALSE  # 简化输出
            )
            
            if (mode_result$converged) {
                theta_robust <- mode_result$theta_mode
                cat(sprintf("  ✓ Weighted Mode: %.4f (SE=%.4f, n=%d IVs)\n", 
                           theta_robust, mode_result$se_mode, length(mode_selected_ivs)))
            } else {
                cat("  ⚠️  Weighted Mode failed\n")
                use_weighted_mode <- FALSE
            }
        } else {
            cat(sprintf("  ⚠️  Insufficient valid IVs (%d)\n", length(mode_selected_ivs)))
            use_weighted_mode <- FALSE
        }
    }
    
    # ========== 策略2: 宽松筛选IV用于先验估计 ==========
    cat("\n[Strategy 2] Liberal IV Selection for Prior Estimation:\n")
    
    # 使用原有的adaptive thresholds逻辑（更宽松）
    thresholds <- determine_adaptive_thresholds(p_X, 2*stats::pnorm(-abs(beta_Y/se_Y)))
    thresh_X_strong <- thresholds$thresh_X_strong
    S_liberal <- which(p_X < thresh_X_strong)
    
    cat(sprintf("  Using p_X < %.2e for priors: %d IVs\n", 
                thresh_X_strong, length(S_liberal)))
    
    # Fallback: 如果Weighted Mode失败，用S_liberal计算IVW作为theta_robust
    if (!use_weighted_mode || is.na(theta_robust)) {
        cat("  → Computing IVW with liberal IVs as fallback theta:\n")
        
        if (length(S_liberal) >= 3) {
            valid_idx <- S_liberal[abs(beta_X[S_liberal]) > 1e-6]
            
            if (length(valid_idx) >= 3) {
                weights <- 1 / (se_Y[valid_idx]^2)
                theta_robust <- sum(weights * beta_Y[valid_idx] * beta_X[valid_idx]) /
                               sum(weights * beta_X[valid_idx]^2)
                
                if (is.finite(theta_robust)) {
                    cat(sprintf("    ✓ IVW: %.4f (n=%d IVs)\n", theta_robust, length(valid_idx)))
                } else {
                    theta_robust <- 0
                    cat("    ⚠️  IVW failed, using theta=0\n")
                }
            } else {
                theta_robust <- 0
            }
        } else {
            theta_robust <- 0
        }
    }
    
    # ========== Sigma_theta_c（基于liberal IVs）==========
    if (length(S_liberal) >= 5) {
        valid_idx <- S_liberal[abs(beta_X[S_liberal]) > 1e-6]
        
        if (length(valid_idx) >= 5) {
            ratios <- abs(beta_Y[valid_idx] / beta_X[valid_idx])
            ratios <- ratios[is.finite(ratios)]
            
            if (length(ratios) >= 5) {
                gamma_95 <- stats::quantile(ratios, 0.95, na.rm = TRUE)
                sigma_theta_c <- gamma_95 / 1.5
                sigma_theta_c <- max(min(sigma_theta_c, 10.0), 0.5)
            } else {
                sigma_theta_c <- 1.5
            }
        } else {
            sigma_theta_c <- 1.5
        }
    } else {
        sigma_theta_c <- 1.5
    }
    
    # 如果theta_robust很大，扩展sigma_theta_c
    if (abs(theta_robust) > sigma_theta_c) {
        sigma_theta_c <- min(abs(theta_robust) * 1.5, 10.0)
        cat(sprintf("  [Adjustment] σ_θ expanded to %.4f (|theta|=%.4f)\n",
                    sigma_theta_c, abs(theta_robust)))
    }
    
    cat(sprintf("  Final sigma_theta_c: %.4f\n", sigma_theta_c))
    
    # ========== 其他先验（使用S_liberal）==========
    S_strong <- S_liberal  # 使用liberal筛选结果
    
    if (length(S_strong) >= 3) {
        sigma_x_sq_hat <- max(stats::var(beta_X[S_strong]) - mean(se_X[S_strong]^2), 1e-6)
        a_x <- 3.0
        b_x <- sigma_x_sq_hat * (a_x - 1)
        b_x <- max(b_x, 0.05)
    } else {
        a_x <- 2.0
        b_x <- 0.1
    }
    cat(sprintf("  Sigma_x_sq Prior: IG(a=%.1f, b=%.4f), E=%.4f\n", 
                a_x, b_x, b_x / (a_x - 1)))
    
    if (length(S_strong) >= 5) {
        weights <- 1 / se_Y[S_strong]^2
        tryCatch({
            fit_egger <- stats::lm(beta_Y[S_strong] ~ beta_X[S_strong], weights = weights)
            coef_egger <- summary(fit_egger)$coefficients
            intercept_p <- coef_egger["(Intercept)", "Pr(>|t|)"]
            
            if (intercept_p < 0.05) {
                sigma_y_sq_hat <- max(summary(fit_egger)$sigma^2, 1e-6)
                a_y <- 2.5
                b_y <- sigma_y_sq_hat * (a_y - 1)
                cat(sprintf("  Sigma_y_sq Prior: UHP DETECTED. IG(a=%.1f, b=%.4f), E=%.4f\n", 
                            a_y, b_y, b_y / (a_y - 1)))
            } else {
                a_y <- 2.5
                b_y <- 0.05 * sigma_theta_c^2 * (a_y - 1)
                cat(sprintf("  Sigma_y_sq Prior: UHP NOT detected. IG(a=%.1f, E=%.4f)\n", 
                            a_y, b_y / (a_y - 1)))
            }
            b_y <- max(b_y, 0.001)
        }, error = function(e) {
            a_y <<- 2.0
            b_y <<- 0.1 * sigma_theta_c^2
            cat(sprintf("  Sigma_y_sq Prior: Egger failed. IG(a=%.1f, b=%.4f)\n", a_y, b_y))
        })
    } else {
        a_y <- 2.0
        b_y <- 0.1 * sigma_theta_c^2
        cat(sprintf("  Sigma_y_sq Prior: Insufficient IVs. IG(a=%.1f, b=%.4f)\n", a_y, b_y))
    }
    
    # ========== 4成分模型：CHP先验（使用S_strong）==========
    if (n_components == 4) {
        cat("\n--- 4-Component CHP Prior Design (liberal IVs) ---\n")
        
        if (length(S_strong) >= 5) {
            ratios <- beta_Y[S_strong] / beta_X[S_strong]
            valid_ratios_idx <- is.finite(ratios) & (abs(beta_X[S_strong]) > 1e-6)
            
            if (sum(valid_ratios_idx) >= 5) {
                weights_ivw <- 1 / se_Y[S_strong[valid_ratios_idx]]^2
                ivw_estimate <- sum(weights_ivw * ratios[valid_ratios_idx]) / sum(weights_ivw)
                
                ratio_sd <- stats::sd(ratios[valid_ratios_idx], na.rm = TRUE)
                q75 <- stats::quantile(abs(ratios[valid_ratios_idx] - ivw_estimate), 0.75, na.rm = TRUE)
                q25 <- stats::quantile(abs(ratios[valid_ratios_idx] - ivw_estimate), 0.25, na.rm = TRUE)
                iqr <- q75 - q25
                
                sd_threshold <- ifelse(iqr / ratio_sd > 0.5, 1.5, 2.0)
                chp_candidates_idx <- which(abs(ratios[valid_ratios_idx] - ivw_estimate) > 
                                           sd_threshold * ratio_sd)
                chp_candidates <- S_strong[valid_ratios_idx][chp_candidates_idx]
                
                if (length(chp_candidates) >= 3) {
                    beta2_est <- stats::median(ratios[chp_candidates], na.rm = TRUE)
                    
                    all_residuals <- beta_Y[S_strong] - beta2_est * beta_X[S_strong]
                    all_se_Y_sq <- se_Y[S_strong]^2
                    
                    mad_residuals <- stats::median(abs(all_residuals - stats::median(all_residuals)))
                    sigma_tilde_c_sq_hat_robust <- (1.4826 * mad_residuals)^2
                    sigma_tilde_c_sq_hat_mom <- max(stats::var(all_residuals) - mean(all_se_Y_sq), 1e-6)
                    
                    sigma_tilde_c_sq_hat <- sqrt(sigma_tilde_c_sq_hat_robust * sigma_tilde_c_sq_hat_mom)
                    
                    a_tilde_c <- 2.5
                    b_tilde_c <- sigma_tilde_c_sq_hat * (a_tilde_c - 1)
                    b_tilde_c <- max(b_tilde_c, 1e-6)
                    
                    cat(sprintf("  [CHP Prior] IG(a=%.2f, b=%.4f), E=%.6f (n=%d IVs)\n", 
                                a_tilde_c, b_tilde_c, b_tilde_c / (a_tilde_c - 1), 
                                length(S_strong)))
                    
                } else {
                    a_tilde_c <- 2.5
                    b_tilde_c <- 0.1 * sigma_theta_c^2 * (a_tilde_c - 1)
                    cat(sprintf("  [CHP Prior] Few candidates, default: IG(a=%.1f, E=%.4f)\n", 
                                a_tilde_c, b_tilde_c / (a_tilde_c - 1)))
                }
                
            } else {
                a_tilde_c <- 2.5
                b_tilde_c <- 0.1 * sigma_theta_c^2 * (a_tilde_c - 1)
            }
        } else {
            a_tilde_c <- 2.5
            b_tilde_c <- 0.1 * sigma_theta_c^2 * (a_tilde_c - 1)
        }
        
        # ========== Summary ==========
        cat("\n[Summary] Dual IV Selection:\n")
        cat(sprintf("  Weighted Mode: %d IVs (p_X < %.2e) → theta_init = %.4f\n",
                   length(mode_selected_ivs), mode_threshold, theta_robust))
        cat(sprintf("  Priors: %d IVs (p_X < %.2e)\n",
                   length(S_liberal), thresh_X_strong))
        
        return(list(
            sigma_theta_c = sigma_theta_c,
            a_x = a_x, b_x = b_x,
            a_y = a_y, b_y = b_y,
            a_tilde_c = a_tilde_c, 
            b_tilde_c = b_tilde_c,
            theta_robust = theta_robust,
            mode_result = mode_result,
            mode_selected_ivs = mode_selected_ivs,      # 严格筛选的IVs
            mode_threshold = mode_threshold,            # 严格阈值
            prior_selected_ivs = S_liberal,             # 宽松筛选的IVs
            prior_threshold = thresh_X_strong           # 宽松阈值
        ))
        
    } else {
        # 3-component model
        cat("\n[Summary] Dual IV Selection:\n")
        cat(sprintf("  Weighted Mode: %d IVs (p_X < %.2e) → theta_init = %.4f\n",
                   length(mode_selected_ivs), mode_threshold, theta_robust))
        cat(sprintf("  Priors: %d IVs (p_X < %.2e)\n",
                   length(S_liberal), thresh_X_strong))
        
        return(list(
            sigma_theta_c = sigma_theta_c,
            a_x = a_x, b_x = b_x,
            a_y = a_y, b_y = b_y,
            theta_robust = theta_robust,
            mode_result = mode_result,
            mode_selected_ivs = mode_selected_ivs,
            mode_threshold = mode_threshold,
            prior_selected_ivs = S_liberal,
            prior_threshold = thresh_X_strong
        ))
    }
}

################################################################################
# PART 8: Check convergence
################################################################################

check_convergence <- function(fit, n_chains = 4, iter = 1000, warmup = 200, n_components = 3) {
  cat("\n===== Convergence Diagnostics =====\n")
  
  summary_fit <- fit$summary()
  
  if (n_components == 3) {
    key_params <- c("theta", "sigma_x_sq", "sigma_y_sq",
                    paste0("pi[", 1:n_components, "]"))
  } else {
    key_params <- c("theta", "beta2", "sigma_x_sq", "sigma_y_sq", "sigma_tilde_c_sq",
                    paste0("pi[", 1:n_components, "]"))
  }
  
  summary_subset <- summary_fit[summary_fit$variable %in% key_params, ]
  
  rhat_vals <- summary_subset$rhat
  names(rhat_vals) <- summary_subset$variable
  n_eff_vals <- summary_subset$ess_bulk
  names(n_eff_vals) <- summary_subset$variable
  
  issues <- character(0)
  
  if (any(is.na(rhat_vals))) {
    issues <- c(issues, "Rhat contains NA")
  }
  if (any(rhat_vals > 1.1, na.rm = TRUE)) {
    bad_params <- names(which(rhat_vals > 1.1))
    issues <- c(issues, sprintf("Rhat > 1.1: %s", paste(bad_params, collapse=", ")))
  }
  
  if (any(n_eff_vals < 100, na.rm = TRUE)) {
    bad_params <- names(which(n_eff_vals < 100))
    issues <- c(issues, sprintf("ESS < 100: %s", paste(bad_params, collapse=", ")))
  }
  
  n_total <- (iter - warmup) * n_chains
  ess_ratios <- n_eff_vals / n_total
  min_ess_ratio <- min(ess_ratios, na.rm = TRUE)
  
  if (any(ess_ratios < 0.1, na.rm = TRUE)) {
    bad_params <- names(which(ess_ratios < 0.1))
    issues <- c(issues, sprintf("ESS/Total < 0.1: %s", paste(bad_params, collapse=", ")))
  }
  
  divergences <- fit$diagnostic_summary()$num_divergent
  total_divergences <- sum(divergences)
  
  if (total_divergences > 0) {
    issues <- c(issues, sprintf("%d divergent transitions", total_divergences))
  }
  
  max_treedepth <- fit$diagnostic_summary()$num_max_treedepth
  if (sum(max_treedepth) > 0) {
    issues <- c(issues, "Maximum tree depth reached")
  }
  
  converged <- length(issues) == 0
  
  cat(sprintf("\nParameter diagnostics:\n"))
  cat(sprintf("  Rhat range: [%.4f, %.4f]\n",
              min(rhat_vals, na.rm = TRUE), max(rhat_vals, na.rm = TRUE)))
  cat(sprintf("  ESS range: [%.0f, %.0f]\n",
              min(n_eff_vals, na.rm = TRUE), max(n_eff_vals, na.rm = TRUE)))
  cat(sprintf("  Min ESS ratio: %.3f\n", min_ess_ratio))
  
  cat(sprintf("\nSampler diagnostics:\n"))
  cat(sprintf("  Divergent transitions: %d\n", total_divergences))
  cat(sprintf("  Max tree depth warnings: %d\n", sum(max_treedepth)))
  
  if (!converged) {
    cat("\n[WARNING] Convergence issues detected:\n")
    for (issue in issues) {
      cat(sprintf("  - %s\n", issue))
    }
  } else {
    cat("\n[OK] All diagnostics passed\n")
  }
  
  return(list(
    converged = converged,
    issues = issues,
    divergences = total_divergences,
    min_ess_ratio = min_ess_ratio,
    max_treedepth = sum(max_treedepth)
  ))
}

################################################################################
# PART 9: Set initial values based on Weighted Mode
################################################################################

get_initial_values <- function(beta_X, beta_Y, se_X, se_Y, priors, 
                               n_chains = 4, n_components = 3) {
  
  valid_idx <- is.finite(beta_X) & is.finite(beta_Y) & se_X > 0 & se_Y > 0
  
  # 优先使用先验中的mode结果
  if (!is.null(priors$mode_result) && priors$mode_result$converged) {
    theta_center <- priors$mode_result$theta_mode
    cat(sprintf("  [Initialization] Using Weighted Mode: %.4f\n", theta_center))
  } else if (!is.null(priors$theta_robust) && is.finite(priors$theta_robust)) {
    theta_center <- priors$theta_robust
    cat(sprintf("  [Initialization] Using robust theta: %.4f\n", theta_center))
  } else {
    # Fallback: IVW
    if (sum(valid_idx) >= 5) {
      weights <- 1 / (se_Y[valid_idx]^2)
      theta_center <- sum(weights * beta_Y[valid_idx] * beta_X[valid_idx]) /
                      sum(weights * beta_X[valid_idx]^2)
      
      if (!is.finite(theta_center)) {
        theta_center <- 0
      } else {
        max_init <- 5 * priors$sigma_theta_c
        theta_center <- max(min(theta_center, max_init), -max_init)
      }
    } else {
      theta_center <- 0
    }
    cat(sprintf("  [Initialization] Using IVW: %.4f\n", theta_center))
  }
  
  init_list <- lapply(1:n_chains, function(chain_id) {
    theta_init <- theta_center + stats::rnorm(1, 0, 0.15 * priors$sigma_theta_c)
    max_theta <- 5 * priors$sigma_theta_c
    theta_init <- max(min(theta_init, max_theta), -max_theta)
    
    if (n_components == 3) {
      list(
        theta_raw = theta_init / priors$sigma_theta_c,
        sigma_x_sq = abs(priors$b_x / priors$a_x * rlnorm(1, 0, 0.2)),
        sigma_y_sq = abs(priors$b_y / priors$a_y * rlnorm(1, 0, 0.2)),
        pi = as.numeric(MCMCpack::rdirichlet(1, priors$alpha_pi))
      )
    } else {
      beta2_init <- theta_center + stats::rnorm(1, 0, 0.5 * priors$sigma_theta_c)
      beta2_init <- max(min(beta2_init, max_theta), -max_theta)
      
      list(
        theta_raw = theta_init / priors$sigma_theta_c,
        beta2_raw = beta2_init / priors$sigma_theta_c,
        sigma_x_sq = abs(priors$b_x / priors$a_x * rlnorm(1, 0, 0.2)),
        sigma_y_sq = abs(priors$b_y / priors$a_y * rlnorm(1, 0, 0.2)),
        sigma_tilde_c_sq = abs(priors$b_tilde_c / priors$a_tilde_c * rlnorm(1, 0, 0.3)),
        pi = as.numeric(MCMCpack::rdirichlet(1, priors$alpha_pi))
      )
    }
  })
  
  return(init_list)
}

################################################################################
# PART 10: Fit single model
################################################################################

fit_single_model <- function(beta_X, beta_Y, se_X, se_Y, rho_e = 0,
                             p_X = NULL, p_Y = NULL,
                             pval_threshold_X = NULL,
                             pval_threshold_Y = NULL,
                             prior_strength = 1.0,
                             n_components = 3,
                             compiled_model = NULL,
                             n_chains = 4, iter = 1000, warmup = 200,
                             use_weighted_mode = TRUE) {
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Fitting %d-component model\n", n_components))
  cat(sprintf("========================================\n"))
  
  J <- length(beta_X)
  
  # 计算先验（包含weighted mode）
  dirichlet_info <- compute_smart_dirichlet_prior(
    beta_X, beta_Y, se_X, se_Y, p_X, p_Y,
    pval_threshold_X, pval_threshold_Y, prior_strength,
    n_components = n_components
  )
  
  other_priors <- compute_empirical_priors(beta_X, beta_Y, se_X, se_Y, p_X, 
                                          n_components, use_weighted_mode)
  priors <- c(other_priors, list(alpha_pi = dirichlet_info$alpha_pi))
  
  cat(sprintf("\n  sigma(theta) = %.4f\n", priors$sigma_theta_c))
  cat(sprintf("  E[sigma^2_x] = %.6f\n", priors$b_x / priors$a_x))
  cat(sprintf("  E[sigma^2_y] = %.6f\n", priors$b_y / priors$a_y))
  if (n_components == 4) {
    cat(sprintf("  E[sigma^2_tilde_c] = %.6f\n", priors$b_tilde_c / priors$a_tilde_c))
  }
  
  # 准备Stan数据
  stan_data <- list(
    J = J, beta_X_hat = beta_X, beta_Y_hat = beta_Y,
    se_X = se_X, se_Y = se_Y, rho_e = rho_e,
    sigma_theta_c = priors$sigma_theta_c,
    a_x = priors$a_x, b_x = priors$b_x,
    a_y = priors$a_y, b_y = priors$b_y,
    alpha_pi = priors$alpha_pi
  )
  
  if (n_components == 4) {
    stan_data$a_tilde_c <- priors$a_tilde_c
    stan_data$b_tilde_c <- priors$b_tilde_c
  }
  
  # 初始值
  init_values <- get_initial_values(beta_X, beta_Y, se_X, se_Y, priors, n_chains, n_components)
  
  # 拟合模型
  cat(sprintf("\nFitting Stan model with cmdstanr...\n"))
  fit <- compiled_model$sample(
    data = stan_data,
    init = init_values,
    chains = n_chains,
    parallel_chains = min(n_chains, parallel::detectCores()),
    iter_warmup = warmup,
    iter_sampling = iter - warmup,
    refresh = 0,
    show_messages = FALSE
  )
  
  # 收敛诊断
  conv_result <- check_convergence(fit, n_chains, iter, warmup, n_components)
  
  # 提取结果
  draws <- fit$draws(format = "df")
  
  theta_samples <- draws$theta
  theta_mean <- mean(theta_samples)
  theta_se <- stats::sd(theta_samples)
  theta_ci <- stats::quantile(theta_samples, c(0.025, 0.975))
  
  z_wald <- theta_mean / theta_se
  p_wald <- 2 * stats::pnorm(-abs(z_wald))
  
  # 提取pi
  pi_cols <- paste0("pi[", 1:n_components, "]")
  pi_samples <- suppressWarnings(as.matrix((draws[, pi_cols])))
  pi_mean <- colMeans(pi_samples)
  pi_sd <- apply(pi_samples, 2, sd)
  
  # 对于4成分模型
  if (n_components == 4) {
    beta2_samples <- draws$beta2
    delta_samples <- draws$delta
    
    beta2_mean <- mean(beta2_samples)
    beta2_ci <- stats::quantile(beta2_samples, c(0.025, 0.975))
    
    delta_mean <- mean(delta_samples)
    delta_ci <- stats::quantile(delta_samples, c(0.025, 0.975))
  } else {
    beta2_samples <- NULL
    delta_samples <- NULL
    beta2_mean <- NA
    beta2_ci <- c(NA, NA)
    delta_mean <- NA
    delta_ci <- c(NA, NA)
  }
  
  # 计算LOO
  log_lik_cols <- paste0("log_lik[", 1:J, "]")
  log_lik_matrix <- suppressWarnings(as.matrix(draws[, log_lik_cols]))
  loo_result <- compute_robust_loo(fit, log_lik_matrix)
  
  # 如果使用了weighted mode，对比结果
  if (!is.null(priors$mode_result) && priors$mode_result$converged) {
    cat("\n=== Comparison: Weighted Mode vs Bayesian ===\n")
    cat(sprintf("  Weighted Mode: %.4f (SE=%.4f)\n", 
                priors$mode_result$theta_mode, priors$mode_result$se_mode))
    cat(sprintf("  Bayesian: %.4f (SE=%.4f)\n", theta_mean, theta_se))
    
    diff <- abs(priors$mode_result$theta_mode - theta_mean)
    if (diff < 2 * theta_se) {
      cat("  ✓ Estimates are consistent\n")
    } else {
      cat("  ⚠️  Substantial difference - check convergence\n")
    }
  }
  
  cat(sprintf("\n%d-component model fitted successfully\n", n_components))
  cat(sprintf("  theta = %.4f (95%% CI: [%.4f, %.4f])\n", theta_mean, theta_ci[1], theta_ci[2]))
  
  result <- list(
    fit = fit,
    priors = priors,
    dirichlet_info = dirichlet_info,
    convergence = conv_result,
    n_components = n_components,
    pi_posterior = list(mean = pi_mean, sd = pi_sd),
    theta = list(
      mean = theta_mean,
      se = theta_se,
      ci = theta_ci,
      z = z_wald,
      p_value = p_wald,
      samples = theta_samples
    ),
    loo_result = loo_result,
    mode_result = priors$mode_result
  )
  
  if (n_components == 4) {
    result$beta2 <- list(
      mean = beta2_mean,
      ci = beta2_ci,
      samples = beta2_samples
    )
    result$delta <- list(
      mean = delta_mean,
      ci = delta_ci,
      samples = delta_samples
    )
  }
  
  return(result)
}

################################################################################
# PART 11: Compare models
################################################################################

compare_models_elpd <- function(result_3comp, result_4comp,
                               significance_level = 0.05) {
  
  cat("\n\n")
  cat("================================================================================\n")
  cat("                          MODEL COMPARISON VIA ELPD                            \n")
  cat("================================================================================\n\n")
  
  loo_3 <- result_3comp$loo_result
  loo_4 <- result_4comp$loo_result
  
  comp <- loo::loo_compare(list(
    "3-comp" = loo_3,
    "4-comp" = loo_4
  ))
  
  cat("LOO comparison table:\n")
  cat("(Models ranked by ELPD, best model first)\n\n")
  print(comp)
  
  best_model_name <- rownames(comp)[1]
  
  cat("\n\n")
  cat("================================================================================\n")
  cat("                        PAIRWISE SIGNIFICANCE TESTS                            \n")
  cat("================================================================================\n\n")
  cat("Using one-sided normal test: Z = ΔELPD / SE(ΔELPD)\n")
  cat(sprintf("Significance level: α = %.3f (Z_crit = %.3f)\n\n",
              significance_level, stats::qnorm(1 - significance_level)))
  
  comparisons <- list()
  
  for (i in 2:nrow(comp)) {
    model_i <- rownames(comp)[i]
    elpd_diff <- comp[i, "elpd_diff"]
    se_diff <- comp[i, "se_diff"]
    
    z_stat <- abs(elpd_diff) / se_diff
    p_value_one_sided <- stats::pnorm(-z_stat)
    
    is_significant <- (z_stat > stats::qnorm(1 - significance_level))
    
    cat(sprintf("Compare: %s vs %s\n", best_model_name, model_i))
    cat(sprintf("  ΔELPD = %.2f (SE = %.2f)\n", elpd_diff, se_diff))
    cat(sprintf("  Z = %.3f, one-sided p = %.4f %s\n",
                z_stat, p_value_one_sided,
                ifelse(is_significant, "***", "")))
    cat(sprintf("  Conclusion: %s\n\n",
                ifelse(is_significant,
                      sprintf("%s is significantly better", best_model_name),
                      "No significant difference")))
    
    comparisons[[model_i]] <- list(
      elpd_diff = elpd_diff,
      se_diff = se_diff,
      z_stat = z_stat,
      p_value = p_value_one_sided,
      is_significant = is_significant
    )
  }
  
  cat("================================================================================\n")
  cat("                            MODEL SELECTION                                    \n")
  cat("================================================================================\n\n")
  
  if (best_model_name == "3-comp") {
    selected_model <- "3-comp"
    selected_result <- result_3comp
    cat("Selected model: 3-component\n")
    cat("Reason: Simpler model with best ELPD\n")
    
  } else if (best_model_name == "4-comp") {
    if (comparisons[["3-comp"]]$is_significant) {
      selected_model <- "4-comp"
      selected_result <- result_4comp
      cat("Selected model: 4-component\n")
      cat("Reason: Significantly better than 3-component model\n")
    } else {
      selected_model <- "3-comp"
      selected_result <- result_3comp
      cat("Selected model: 3-component (Occam's razor)\n")
      cat("Reason: No significant improvement from 4-component\n")
    }
  }
  
  cat("\n")
  cat("================================================================================\n")
  cat("                         FINAL INFERENCE RESULTS                               \n")
  cat("================================================================================\n\n")
  
  cat(sprintf("Based on %s:\n\n", selected_model))
  cat(sprintf("Causal effect (theta):\n"))
  cat(sprintf("  Estimate: %.4f (SE = %.4f)\n",
              selected_result$theta$mean, selected_result$theta$se))
  cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
              selected_result$theta$ci[1], selected_result$theta$ci[2]))
  cat(sprintf("  P-value: %.4e %s\n",
              selected_result$theta$p_value,
              ifelse(selected_result$theta$p_value < 0.05, "***", "")))
  
  if (selected_model == "4-comp") {
    cat(sprintf("\nConfounding bias (delta = beta2 - theta):\n"))
    cat(sprintf("  Estimate: %.4f\n", selected_result$delta$mean))
    cat(sprintf("  95%% CI: [%.4f, %.4f]\n",
                selected_result$delta$ci[1], selected_result$delta$ci[2]))
    
    if (selected_result$delta$ci[1] > 0) {
      cat("  Interpretation: Significant POSITIVE confounding bias\n")
    } else if (selected_result$delta$ci[2] < 0) {
      cat("  Interpretation: Significant NEGATIVE confounding bias\n")
    } else {
      cat("  Interpretation: No significant confounding bias\n")
    }
  }
  
  cat(sprintf("\nComponent weights (posterior mean):\n"))
  n_comp <- selected_result$n_components
  if (n_comp == 3) {
    comp_labels <- c("Valid", "Valid+UHP", "NULL/UHP")
  } else {
    comp_labels <- c("Valid", "Valid+UHP", "Valid+UHP+CHP", "NULL")
  }
  
  for (i in 1:n_comp) {
    cat(sprintf("  π[%d] (%15s): %.3f (SD = %.3f)\n",
                i, comp_labels[i],
                selected_result$pi_posterior$mean[i],
                selected_result$pi_posterior$sd[i]))
  }
  
  cat("\n")
  cat("================================================================================\n\n")
  
  return(list(
    comparison_table = comp,
    pairwise_tests = comparisons,
    selected_model = selected_model,
    selected_result = selected_result,
    best_elpd_model = best_model_name
  ))
}

################################################################################
# PART 12: Main function
################################################################################

run_mr_model_comparison <- function(beta_X, beta_Y, se_X, se_Y, rho_e = 0,
                                   p_X = NULL, p_Y = NULL,
                                   pval_threshold_X = NULL,
                                   pval_threshold_Y = NULL,
                                   prior_strength = 1.0,
                                   n_chains = 4, iter = 1000, warmup = 200,
                                   significance_level = 0.05,
                                   use_parallel = TRUE,
                                   use_weighted_mode = TRUE,
                                   compiled_models = NULL) {
  
  cat("\n")
  cat("################################################################################\n")
  cat("#    MR MODEL COMPARISON: 3 vs 4 COMPONENTS (with Weighted Mode Init)         #\n")
  cat("################################################################################\n\n")
  
  J <- length(beta_X)
  stopifnot("Need at least 5 SNPs" = J >= 5)
  stopifnot("All se must be positive" = all(se_X > 0) && all(se_Y > 0))
  stopifnot("rho_e must be in [-1, 1]" = abs(rho_e) <= 1)
  
  cat(sprintf("Analysis configuration:\n"))
  cat(sprintf("  Number of SNPs: %d\n", J))
  cat(sprintf("  Sample overlap: ρ = %.3f\n", rho_e))
  cat(sprintf("  MCMC: %d chains, %d iterations (%d warmup)\n", n_chains, iter, warmup))
  cat(sprintf("  Weighted Mode initialization: %s\n", ifelse(use_weighted_mode, "YES", "NO")))
  
  # 编译模型
  if (is.null(compiled_models)) {
    cat("\nCompiling Stan models...\n")
    compiled_models <- compile_stan_models()
  } else {
    cat("\nUsing pre-compiled Stan models\n")
  }
  
  # 并行配置
  total_cores <- parallel::detectCores()
  cat(sprintf("  Available cores: %d\n", total_cores))
  
  if (use_parallel && total_cores >= 4) {
    if (.Platform$OS.type == "windows") {
      future::plan(future::multisession, workers = 2)
      cat("  Parallel mode: ENABLED (future::multisession)\n")
    } else {
      future::plan(future::multicore, workers = 2)
      cat("  Parallel mode: ENABLED (future::multicore)\n")
    }
    
    if (total_cores >= 8) {
      chains_per_model <- n_chains
    } else if (total_cores >= 4) {
      chains_per_model <- max(2, floor(total_cores / 2))
    } else {
      chains_per_model <- 2
    }
    
    if (chains_per_model != n_chains) {
      cat(sprintf("  Adjusted chains: %d -> %d per model\n", n_chains, chains_per_model))
    }
    
  } else {
    future::plan(future::sequential)
    chains_per_model <- n_chains
    cat("  Parallel mode: DISABLED (future::sequential)\n")
  }
  
  cat(sprintf("  Significance level: α = %.3f\n\n", significance_level))
  
  # 拟合模型
  cat("===== Fitting Models =====\n")
  start_time <- Sys.time()
  
  if (use_parallel && !inherits(future::plan(), "future::sequential")) {
    cat("Running 3-comp and 4-comp models in PARALLEL...\n\n")
    
    results_list <- future.apply::future_lapply(
      X = c(3, 4),
      FUN = function(n_comp) {
        compiled_model <- if (n_comp == 3) {
          compiled_models$model_3comp
        } else {
          compiled_models$model_4comp
        }
        
        fit_single_model(
          beta_X = beta_X,
          beta_Y = beta_Y,
          se_X = se_X,
          se_Y = se_Y,
          rho_e = rho_e,
          p_X = p_X,
          p_Y = p_Y,
          pval_threshold_X = pval_threshold_X,
          pval_threshold_Y = pval_threshold_Y,
          prior_strength = prior_strength,
          n_components = n_comp,
          compiled_model = compiled_model,
          n_chains = chains_per_model,
          iter = iter,
          warmup = warmup,
          use_weighted_mode = use_weighted_mode
        )
      },
      future.seed = TRUE
    )
    
    result_3comp <- results_list[[1]]
    result_4comp <- results_list[[2]]
    
  } else {
    cat("Running models SEQUENTIALLY...\n\n")
    
    result_3comp <- fit_single_model(
      beta_X, beta_Y, se_X, se_Y, rho_e,
      p_X, p_Y, pval_threshold_X, pval_threshold_Y,
      prior_strength, n_components = 3,
      compiled_model = compiled_models$model_3comp,
      chains_per_model, iter, warmup, use_weighted_mode
    )
    
    result_4comp <- fit_single_model(
      beta_X, beta_Y, se_X, se_Y, rho_e,
      p_X, p_Y, pval_threshold_X, pval_threshold_Y,
      prior_strength, n_components = 4,
      compiled_model = compiled_models$model_4comp,
      chains_per_model, iter, warmup, use_weighted_mode
    )
  }
  
  elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  cat("\n===== Model Fitting Complete =====\n")
  cat(sprintf("Total elapsed time: %.2f minutes\n", elapsed_time))
  
  # 重置为future::sequential
  future::plan(future::sequential)
  
  # 模型比较
  comparison_results <- compare_models_elpd(
    result_3comp, result_4comp,
    significance_level
  )
  
  return(list(
    result_3comp = result_3comp,
    result_4comp = result_4comp,
    comparison = comparison_results,
    selected_model = comparison_results$selected_model,
    final_result = comparison_results$selected_result,
    elapsed_time_minutes = elapsed_time,
    parallel_used = use_parallel && (total_cores >= 4),
    compiled_models = compiled_models,
    used_weighted_mode = use_weighted_mode
  ))
}

################################################################################
# PART 13: Visualization
################################################################################

plot_model_comparison <- function(comparison_results) {
  
  comp_table <- comparison_results$comparison_table
  
  models <- rownames(comp_table)
  elpd <- comp_table[, "elpd_diff"] + comp_table[1, "elpd_loo"]
  se <- comp_table[, "se_diff"]
  
  se[1] <- comparison_results$selected_result$loo_result$estimates["elpd_loo", "SE"]
  
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))
  
  x_pos <- 1:length(models)
  plot(x_pos, elpd, ylim = range(c(elpd - 2*se, elpd + 2*se)),
       xlab = "", ylab = "ELPD", xaxt = "n",
       main = "Model Comparison via LOO-CV",
       pch = 19, cex = 1.5, col = "steelblue")
  
  arrows(x_pos, elpd - se, x_pos, elpd + se,
         angle = 90, code = 3, length = 0.1, col = "steelblue", lwd = 2)
  
  axis(1, at = x_pos, labels = models, las = 1)
  
  abline(h = elpd[1], lty = 2, col = "red")
  
  best_model <- comparison_results$selected_model
  legend("bottomright",
         legend = c("Best model", sprintf("Selected: %s", best_model)),
         col = c("red", "black"), lty = c(2, NA), pch = c(NA, NA),
         bty = "n", cex = 0.9)
  
  grid()
}

plot_theta_comparison <- function(result_3comp, result_4comp) {
  
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  
  hist(result_3comp$theta$samples, breaks = 50, col = "lightblue", border = "white",
       main = "3-component", xlab = "theta", freq = FALSE, xlim = range(c(
         result_3comp$theta$samples,
         result_4comp$theta$samples
       )))
  abline(v = result_3comp$theta$mean, col = "red", lwd = 2)
  abline(v = 0, col = "black", lty = 2)
  
  hist(result_4comp$theta$samples, breaks = 50, col = "lightgreen", border = "white",
       main = "4-component", xlab = "theta", freq = FALSE, xlim = range(c(
         result_3comp$theta$samples,
         result_4comp$theta$samples
       )))
  abline(v = result_4comp$theta$mean, col = "red", lwd = 2)
  abline(v = 0, col = "black", lty = 2)
  
  par(mfrow = c(1, 1))
}

################################################################################
# PART 14: Simulation 
################################################################################

simulate_mr_data <- function(J = 200, theta_true = 0, delta_true = 0.2,
                             prop_valid = 0.3, prop_valid_uhp = 0.1, 
                             prop_chp = 0.2, prop_chp_uhp = 0.1,
                             prop_uhp = 0.1, prop_null = 0.2) {
  
  props <- c(prop_valid, prop_valid_uhp, prop_chp, prop_chp_uhp, prop_uhp, prop_null)
  props <- props / sum(props)
  
  components <- sample(1:6, J, replace = TRUE, prob = props)
  
  beta2_true <- theta_true + delta_true
  
  u_x <- stats::rnorm(J, 0, 0.05)
  u_y_uhp <- stats::rnorm(J, 0, 0.02)
  tilde_c <- stats::rnorm(J, 0, 0.01)
  
  beta_X_true <- numeric(J)
  beta_Y_true <- numeric(J)
  
  for (j in 1:J) {
    if (components[j] == 1) {
      beta_X_true[j] <- u_x[j]
      beta_Y_true[j] <- theta_true * u_x[j]
    } else if (components[j] == 2) {
      beta_X_true[j] <- u_x[j]
      beta_Y_true[j] <- theta_true * u_x[j] + u_y_uhp[j]
    } else if (components[j] == 3) {
      beta_X_true[j] <- u_x[j]
      beta_Y_true[j] <- beta2_true * u_x[j] + tilde_c[j]
    } else if (components[j] == 4) {
      beta_X_true[j] <- u_x[j]
      beta_Y_true[j] <- beta2_true * u_x[j] + tilde_c[j] + u_y_uhp[j]
    } else if (components[j] == 5) {
      beta_X_true[j] <- 0
      beta_Y_true[j] <- u_y_uhp[j]
    } else {
      beta_X_true[j] <- 0
      beta_Y_true[j] <- 0
    }
  }
  
  se_X <- runif(J, 0.01, 0.03)
  se_Y <- runif(J, 0.01, 0.03)
  
  beta_X_obs <- beta_X_true + stats::rnorm(J, 0, se_X)
  beta_Y_obs <- beta_Y_true + stats::rnorm(J, 0, se_Y)
  
  return(list(
    beta_X = beta_X_obs,
    beta_Y = beta_Y_obs,
    se_X = se_X,
    se_Y = se_Y,
    components_true = components,
    theta_true = theta_true,
    beta2_true = beta2_true,
    delta_true = delta_true
  ))
}

################################################################################
# Finished
################################################################################

cat("\n>>> MR Model Comparison with Weighted Mode (FINAL VERSION) Loaded\n")
cat(">>> Key improvements:\n")
cat("    1. ✓ Removed MR-Mix spike-detection (no 2-comp vs 4-comp conflict)\n")
cat("    2. ✓ Integrated Weighted Mode (Hartwig 2017) for robust initialization\n")
cat("    3. ✓ All other features preserved\n\n")
cat(">>> Usage:\n")
cat("    1. compiled_models <- compile_stan_models()\n")
cat("    2. results <- run_mr_model_comparison(\n")
cat("           beta_X, beta_Y, se_X, se_Y,\n")
cat("           compiled_models = compiled_models,\n")
cat("           use_weighted_mode = TRUE  # Default\n")
cat("       )\n")
cat("    3. plot_model_comparison(results$comparison)\n\n")
cat(">>> Weighted Mode is now the default robust initialization method!\n\n")
