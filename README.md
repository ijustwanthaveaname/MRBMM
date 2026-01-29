# BMM: Bayesian Mixture Model for Mendelian Randomization

[![R-CMD-check](https://github.com/your-username/BMM/workflows/R-CMD-check/badge.svg)](https://github.com/your-username/BMM/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

BMM implements a robust Bayesian mixture model for Mendelian Randomization that explicitly models horizontal pleiotropy with automatic model selection.

### Key Features

✅ **Explicit pleiotropy modeling** - Handles UHP and CHP with 3/4-component models  
✅ **Robust initialization** - Weighted Mode estimation (Hartwig et al. 2017)  
✅ **Automatic model selection** - LOO-CV based comparison  
✅ **Comprehensive harmonization** - Local plink-based proxy SNP search  
✅ **Sample overlap correction** - Support for overlapping samples  
✅ **Parallel processing** - Multi-core MCMC sampling  

## Installation

### Prerequisites

```r
# 1. Install cmdstanr
install.packages("cmdstanr", 
                repos = c("https://mc-stan.org/r-packages/", 
                         getOption("repos")))

# 2. Install CmdStan
cmdstanr::install_cmdstan()

# 3. Install plink (system requirement)
# Ubuntu/Debian: sudo apt-get install plink
# macOS: brew install plink
# Or download from: https://www.cog-genomics.org/plink/
```

### Install BMM

```r
# From GitHub
devtools::install_github("your-username/BMM")

# Or from source
R CMD INSTALL BMM_1.0.0.tar.gz
```

## Quick Start

```r
library(BMM)

# ⚠️ CRITICAL: Set environment variable
Sys.setenv(TBB_CXX_TYPE = "gcc")

# Paths to plink and reference panel
PLINK_BIN <- "/usr/bin/plink"
REF_PANEL <- "/data/1kg_eur/EUR"  # Without .bed/.bim/.fam extension

# Read GWAS summary statistics
exposure <- read_gwas_data(
  file_path = "exposure.txt",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  ea_col = "EA",
  oa_col = "OA"
)

outcome <- read_gwas_data(
  file_path = "outcome.txt",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  pval_col = "P",
  ea_col = "EA",
  oa_col = "OA"
)

# Harmonize with local plink
harmonized <- harmonize_for_mr(
  exp_gwas = exposure,
  otc_gwas = outcome,
  exp_name = "BMI",
  otc_name = "CAD",
  plink_bin = PLINK_BIN,
  ref_panel = REF_PANEL
)

# Extract data for BMM
bmm_data <- extract_bmm_data(harmonized)

# Compile Stan models (required each session, takes 1-2 min)
compiled <- compile_stan_models()

# Run BMM analysis
results <- run_mr_model_comparison(
  beta_X = bmm_data$beta_X,
  beta_Y = bmm_data$beta_Y,
  se_X = bmm_data$se_X,
  se_Y = bmm_data$se_Y,
  rho_e = 0,  # No sample overlap
  compiled_models = compiled
)

# View results
summary(results)
```

## Output

```
══════════════════════════════════════════════════════════
 BMM Analysis Results
══════════════════════════════════════════════════════════

Selected model: 3-comp
Number of SNPs used: 127

─── Causal Effect (theta) ───
  Estimate: 0.1234
  Std Error: 0.0456
  95% CI: [0.0340, 0.2128]
  P-value: 6.82e-03 ***
  → Significant POSITIVE causal effect

─── Component Weights (Posterior Mean) ───
  π[1] (          Valid): 0.687 (SD = 0.124)
  π[2] (     Valid+UHP): 0.245 (SD = 0.108)
  π[3] (       NULL/UHP): 0.068 (SD = 0.043)

─── Convergence Diagnostics ───
  ✓ All chains converged successfully
    Min ESS ratio: 0.456

══════════════════════════════════════════════════════════
```

## Documentation

- [Complete Workflow](examples/complete_workflow.R) - End-to-end analysis
- [Function Reference](https://your-username.github.io/BMM/reference/) - All functions
- User Guide (coming soon)

## Method

BMM uses a Bayesian mixture model to classify SNPs into:

**3-Component Model:**
1. **Valid IVs** - No pleiotropy (βY = θ·βX)
2. **Valid + UHP** - Uncorrelated horizontal pleiotropy
3. **NULL/UHP** - Null or invalid IVs

**4-Component Model:**
1. **Valid IVs** - No pleiotropy
2. **Valid + UHP** - Uncorrelated horizontal pleiotropy
3. **Valid + UHP + CHP** - Correlated horizontal pleiotropy  
4. **NULL** - Null IVs

Model selection via LOO-CV. Weighted Mode provides robust initialization.

## Citation

```
[Your citation will go here]

For Weighted Mode:
Hartwig FP, Davey Smith G, Bowden J. Robust inference in summary data 
Mendelian randomization via the zero modal pleiotropy assumption. 
Int J Epidemiol. 2017;46(6):1985-1998.
```

## License

MIT License - see LICENSE file

## Contact

- Issues: https://github.com/your-username/BMM/issues
- Email: your.email@university.edu
