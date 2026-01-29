################################################################################
# Complete BMM Workflow Example
# From GWAS data to final results
################################################################################

library(BMM)
library(tidyverse)

# =============================================================================
# ⭐ CRITICAL: Set environment and paths
# =============================================================================

# TBB environment (REQUIRED before Stan)
Sys.setenv(TBB_CXX_TYPE = "gcc")

# Plink and reference panel paths (UPDATE THESE!)
PLINK_BIN <- "/usr/bin/plink"
REF_PANEL <- "/data/1kg_eur/EUR"

# Verify paths
if (!file.exists(PLINK_BIN)) {
  stop("Plink not found: ", PLINK_BIN, "\nPlease update PLINK_BIN!")
}
if (!file.exists(paste0(REF_PANEL, ".bed"))) {
  stop("Reference panel not found: ", REF_PANEL, ".bed\nPlease update REF_PANEL!")
}

cat("✓ Plink:", PLINK_BIN, "\n")
cat("✓ Reference panel:", REF_PANEL, "\n\n")

# =============================================================================
# STEP 1: READ GWAS DATA
# =============================================================================

cat("═══════════════════════════════════════════\n")
cat("STEP 1: Reading GWAS Data\n")
cat("═══════════════════════════════════════════\n\n")

# Example: GWAS with beta/SE
exposure <- read_gwas_data(
  file_path = "data/exposure.txt",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  ea_col = "EA",
  oa_col = "OA",
  eaf_col = "EAF"
)

# Example: GWAS with OR/CI
outcome <- read_gwas_data(
  file_path = "data/outcome.txt",
  snp_col = "SNP",
  pval_col = "P_VALUE",
  ea_col = "RISK_ALLELE",
  oa_col = "OTHER_ALLELE",
  or_col = "OR",
  or_lower_col = "OR_95L",
  or_upper_col = "OR_95U"
)

cat(sprintf("Exposure: %d variants\n", nrow(exposure)))
cat(sprintf("Outcome: %d variants\n\n", nrow(outcome)))

# =============================================================================
# STEP 2: HARMONIZATION (with local plink)
# =============================================================================

cat("═══════════════════════════════════════════\n")
cat("STEP 2: Harmonization\n")
cat("═══════════════════════════════════════════\n\n")

harmonized <- harmonize_for_mr(
  exp_gwas = exposure,
  otc_gwas = outcome,
  exp_name = "BMI",
  otc_name = "CAD",
  
  # Local plink paths
  plink_bin = PLINK_BIN,
  ref_panel = REF_PANEL,
  
  # Parameters
  pval_threshold = 5e-8,
  clump_kb = 10000,
  clump_r2 = 0.001,
  use_proxy = TRUE,
  proxy_r2 = 0.8,
  verbose = TRUE
)

cat(sprintf("\n✓ Harmonized: %d SNPs\n\n", sum(harmonized$mr_keep)))

# =============================================================================
# STEP 3: EXTRACT BMM DATA
# =============================================================================

cat("═══════════════════════════════════════════\n")
cat("STEP 3: Extract Data for BMM\n")
cat("═══════════════════════════════════════════\n\n")

bmm_data <- extract_bmm_data(harmonized)

cat(sprintf("✓ %d SNPs ready for analysis\n\n", bmm_data$n_snps))

# =============================================================================
# STEP 4: COMPILE MODELS
# =============================================================================

cat("═══════════════════════════════════════════\n")
cat("STEP 4: Compile Stan Models\n")
cat("═══════════════════════════════════════════\n\n")

# This step is REQUIRED and takes 1-2 minutes
compiled_models <- compile_bmm_models(
  force_recompile = FALSE,
  set_tbb_env = TRUE
)

# =============================================================================
# STEP 5: RUN BMM
# =============================================================================

cat("═══════════════════════════════════════════\n")
cat("STEP 5: Run BMM Analysis\n")
cat("═══════════════════════════════════════════\n\n")

results <- run_bmm_analysis(
  beta_X = bmm_data$beta_X,
  beta_Y = bmm_data$beta_Y,
  se_X = bmm_data$se_X,
  se_Y = bmm_data$se_Y,
  rho_e = 0,  # No sample overlap
  
  # MCMC settings
  n_chains = 4,
  iter = 1000,
  warmup = 200,
  
  # Options
  use_weighted_mode = TRUE,
  use_parallel = TRUE,
  
  # Pre-compiled models
  compiled_models = compiled_models
)

# =============================================================================
# STEP 6: VIEW RESULTS
# =============================================================================

cat("═══════════════════════════════════════════\n")
cat("STEP 6: Results\n")
cat("═══════════════════════════════════════════\n\n")

summary(results)

# =============================================================================
# STEP 7: SAVE
# =============================================================================

cat("═══════════════════════════════════════════\n")
cat("STEP 7: Save Results\n")
cat("═══════════════════════════════════════════\n\n")

dir.create("output", showWarnings = FALSE)

saveRDS(results, "output/bmm_results.rds")
saveRDS(harmonized, "output/harmonized_data.rds")

cat("✓ Saved: output/bmm_results.rds\n")
cat("✓ Saved: output/harmonized_data.rds\n\n")

cat("✓ Analysis complete!\n\n")
