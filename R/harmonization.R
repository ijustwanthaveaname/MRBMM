################################################################################
# GWAS Data Harmonization - Matching Real Data Evaluation Method
# Uses local plink + reference panel for proxy SNP search
################################################################################
#' Find Proxy SNPs Using Local Plink and Reference Panel for Data Harmonization
#'
#' Exactly matches the get_proxy function from real data evaluation
#' Uses local plink to find LD proxies for missing SNPs
#'
#' @importFrom magrittr %>%
#' @param exp_dat Exposure data (TwoSampleMR formatted)
#' @param otc_dat Outcome data (TwoSampleMR formatted)
#' @param r2 Minimum LD r² threshold (default: 0.8)
#' @param bin_plink Path to plink executable
#' @param ref_panel Path to reference panel (without .bed/.bim/.fam extension)
#' @param outfile Output file prefix for temporary files (default: "proxy")
#'
#' @return Outcome data frame with proxy SNPs added
#'
#' @keywords internal
get_proxy <- function(exp_dat, otc_dat, r2 = 0.8,
                      bin_plink, ref_panel, outfile = "proxy") {
  bin_plink <- normalizePath(bin_plink, mustWork = TRUE)
  ref_panel <- path.expand(ref_panel)
  cat("  Searching for proxy SNPs...\n")
  exp.snp <- exp_dat$SNP
  otc.snp <- otc_dat$SNP

  # Initialize proxy columns
  otc_dat <- otc_dat %>% dplyr::mutate(
    target_snp.outcome = NA,
    proxy_snp.outcome = NA,
    target_a1.outcome = NA,
    target_a2.outcome = NA,
    proxy_a1.outcome = NA,
    proxy_a2.outcome = NA,
    proxy.outcome = FALSE
  )

  # Find missing SNPs
  miss.expsnps <- setdiff(exp.snp, intersect(exp.snp, otc.snp))
  remain.otcsnps <- setdiff(otc.snp, intersect(exp.snp, otc.snp))

  cat(sprintf("    Found %d SNPs missing in outcome\n", length(miss.expsnps)))

  if (length(miss.expsnps) == 0) {
    cat("    No proxy search needed\n")
    return(otc_dat %>% dplyr::filter(SNP %in% exp.snp))
  }

  proxy_count <- 0

  # Search for proxy for each missing SNP
  for (miss.snp in miss.expsnps) {
    # Run plink to find LD proxies
    outfilepath <- shQuote(outfile)
    system(
      glue::glue(
        "{bin_plink} --bfile {ref_panel} --r2 --ld-snp {miss.snp} ",
        "--ld-window 1000000 --ld-window-kb 1000 --ld-window-r2 {r2} ",
        "--out {outfilepath}_{miss.snp}_proxy"
      ),
      ignore.stdout = TRUE, ignore.stderr = TRUE
    )

    # Read LD results
    ld_file <- glue::glue("{outfile}_{miss.snp}_proxy.ld")
    if (!file.exists(ld_file)) {
      next
    }

    proxies.df <- withCallingHandlers(readr::read_table(ld_file, show_col_types = FALSE), warning = function(w) {
      msg <- w$message
      # 忽略 X8 列名填充警告
      if (grepl("Missing column names filled in: 'X8'", msg, fixed = TRUE)) {
        invokeRestart("muffleWarning")
        }
      }
    )

    # Select best proxy (highest R2) that's in outcome data
    top_proxies.df <- proxies.df %>%
      dplyr::filter(SNP_A != SNP_B, SNP_B %in% remain.otcsnps) %>%
      dplyr::arrange(-R2) %>%
      head(1) %>%
      dplyr::rename(RS_Number = SNP_B)

    # Clean up temporary files
    file.remove(ld_file)
    log_file <- glue::glue("{outfile}_{miss.snp}_proxy.log")
    if (file.exists(log_file)) file.remove(log_file)

    if (nrow(top_proxies.df) == 0) {
      next
    }

    # Found a proxy!
    proxy_count <- proxy_count + 1
    rsid <- top_proxies.df$RS_Number

    # Get alleles from exposure
    effect_allele.exposure <- (exp_dat %>% dplyr::filter(SNP == miss.snp) %>% head(1))$effect_allele.exposure
    other_allele.exposure <- (exp_dat %>% dplyr::filter(SNP == miss.snp) %>% head(1))$other_allele.exposure

    # Get allele information from plink
    allele_info <- system(
      glue::glue(
        "{bin_plink} --bfile {ref_panel} --ld {miss.snp} {rsid} | ",
        'grep "In phase" | awk \'{{print substr($NF,1,1)"\\n"substr($NF,4,1)"\\n"',
        'substr($NF,2,1)"\\n"substr($NF,5,1)}}\''
      ),
      intern = TRUE
    )

    if (length(allele_info) != 4) {
      next
    }

    target.allele1 <- allele_info[1]
    target.allele2 <- allele_info[2]
    proxy.allele1 <- allele_info[3]
    proxy.allele2 <- allele_info[4]

    # Check allele consistency
    if (!identical(
      sort(c(target.allele1, target.allele2)),
      sort(c(effect_allele.exposure, other_allele.exposure))
    )) {
      next
    }

    # Update outcome data with proxy information
    otc_dat[otc_dat$SNP == rsid, ]$target_snp.outcome <- miss.snp
    otc_dat[otc_dat$SNP == rsid, ]$proxy_snp.outcome <- rsid
    otc_dat[otc_dat$SNP == rsid, ]$target_a1.outcome <-
      ifelse(otc_dat[otc_dat$SNP == rsid, ]$effect_allele.outcome == proxy.allele1,
        target.allele1, target.allele2
      )
    otc_dat[otc_dat$SNP == rsid, ]$target_a2.outcome <-
      ifelse(otc_dat[otc_dat$SNP == rsid, ]$other_allele.outcome == proxy.allele2,
        target.allele2, target.allele1
      )
    otc_dat[otc_dat$SNP == rsid, ]$proxy_a1.outcome <-
      ifelse(otc_dat[otc_dat$SNP == rsid, ]$effect_allele.outcome == proxy.allele1,
        proxy.allele1, proxy.allele2
      )
    otc_dat[otc_dat$SNP == rsid, ]$proxy_a2.outcome <-
      ifelse(otc_dat[otc_dat$SNP == rsid, ]$other_allele.outcome == proxy.allele2,
        proxy.allele2, proxy.allele1
      )
    otc_dat[otc_dat$SNP == rsid, ]$proxy.outcome <- TRUE
    otc_dat[otc_dat$SNP == rsid, ]$effect_allele.outcome <-
      ifelse(otc_dat[otc_dat$SNP == rsid, ]$effect_allele.outcome == proxy.allele1,
        target.allele1, target.allele2
      )
    otc_dat[otc_dat$SNP == rsid, ]$other_allele.outcome <-
      ifelse(otc_dat[otc_dat$SNP == rsid, ]$other_allele.outcome == proxy.allele2,
        target.allele2, target.allele1
      )
    otc_dat[otc_dat$SNP == rsid, ]$SNP <- miss.snp
  }

  cat(sprintf("    Successfully found %d proxy SNPs\n", proxy_count))

  return(otc_dat %>% dplyr::filter(SNP %in% exp.snp))
}

#' Read GWAS Summary Statistics
#'
#' Read and standardize GWAS data matching real data evaluation format
#'
#' @param file_path Path to GWAS file
#' @param snp_col Column name for SNP ID
#' @param beta_col Column name for beta (can be NA if using OR)
#' @param se_col Column name for SE (can be NA if using OR CI)
#' @param pval_col Column name for p-value
#' @param ea_col Column name for effect allele
#' @param oa_col Column name for other allele
#' @param eaf_col Column name for EAF (can be NA)
#' @param chr_col Column name for chromosome (can be NA)
#' @param pos_col Column name for position (can be NA)
#' @param or_col Column name for OR (if beta_col is NA)
#' @param or_lower_col Column name for OR lower CI (if beta_col is NA)
#' @param or_upper_col Column name for OR upper CI (if beta_col is NA)
#'
#' @return Standardized data frame with columns: SNP, EA, OA, BETA, SE, P, EAF, CHR, POS
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # GWAS with beta and SE
#' exposure <- read_gwas_data(
#'   "exposure.txt",
#'   snp_col = "rsid",
#'   beta_col = "beta",
#'   se_col = "se",
#'   pval_col = "pval",
#'   ea_col = "effect_allele",
#'   oa_col = "other_allele",
#'   eaf_col = "eaf"
#' )
#'
#' # GWAS with OR and CI
#' outcome <- read_gwas_data(
#'   "outcome.txt",
#'   snp_col = "SNP",
#'   pval_col = "P_VALUE",
#'   ea_col = "RISK_ALLELE",
#'   oa_col = "OTHER_ALLELE",
#'   or_col = "OR",
#'   or_lower_col = "OR_95L",
#'   or_upper_col = "OR_95U"
#' )
#' }
read_gwas_data <- function(file_path,
                           snp_col, beta_col = NA, se_col = NA, pval_col,
                           ea_col, oa_col, eaf_col = NA,
                           chr_col = NA, pos_col = NA,
                           or_col = NA, or_lower_col = NA, or_upper_col = NA) {
  cat(sprintf("Reading GWAS from: %s\n", file_path))

  data <- data.table::fread(file_path)
  cat(sprintf("  Raw data: %d variants\n", nrow(data)))

  # Build standardized data frame
  std_data <- data.frame(
    SNP = data[[snp_col]],
    EA = data[[ea_col]],
    OA = data[[oa_col]],
    stringsAsFactors = FALSE
  )

  # Handle beta/SE or OR/CI
  if (!is.na(beta_col)) {
    std_data$BETA <- as.numeric(data[[beta_col]])
    std_data$SE <- as.numeric(data[[se_col]])
  } else if (!is.na(or_col)) {
    or_val <- as.numeric(data[[or_col]])
    or_lower <- as.numeric(data[[or_lower_col]])
    or_upper <- as.numeric(data[[or_upper_col]])
    std_data$BETA <- log(or_val)
    if (!is.na(se_col)) {
      std_data$SE <- as.numeric(data[[se_col]])
    } else {
      std_data$SE <- (log(or_upper) - log(or_lower)) / (2 * 1.96)
    }
  } else {
    stop("Must provide either beta_col/se_col or or_col/or_lower_col/or_upper_col")
  }

  # P-value
  std_data$P <- as.numeric(data[[pval_col]])

  # EAF (optional)
  if (!is.na(eaf_col)) {
    std_data$EAF <- as.numeric(data[[eaf_col]])
  } else {
    std_data$EAF <- NA
  }

  # CHR and POS (optional)
  if (!is.na(chr_col)) {
    std_data$CHR <- data[[chr_col]]
    std_data$POS <- as.numeric(data[[pos_col]])
  } else {
    std_data$CHR <- NA
    std_data$POS <- NA
  }

  # Quality control
  std_data <- std_data %>%
    dplyr::filter(!is.na(SNP), !is.na(BETA), !is.na(SE), !is.na(P)) %>%
    dplyr::filter(SNP != ".", SNP != "", nchar(SNP) > 0) %>%
    dplyr::filter(is.finite(BETA), is.finite(SE), is.finite(P)) %>%
    dplyr::filter(SE > 0, P > 0, P <= 1) %>%
    dplyr::filter(stringr::str_detect(SNP, "^rs\\d+"))

  cat(sprintf("  Clean data: %d variants\n", nrow(std_data)))

  return(std_data)
}

#' Complete GWAS Harmonization Workflow
#'
#' Performs complete harmonization matching real data evaluation:
#' 1. Read GWAS data
#' 2. Filter exposure SNPs by p-value
#' 3. Format for TwoSampleMR
#' 4. LD clumping with local reference panel
#' 5. Proxy SNP search with local plink
#' 6. Harmonization
#'
#' @param exp_gwas Exposure GWAS data frame (from read_gwas_data)
#' @param otc_gwas Outcome GWAS data frame (from read_gwas_data)
#' @param exp_name Exposure trait name
#' @param otc_name Outcome trait name
#' @param plink_bin Path to plink executable (e.g., "/usr/bin/plink")
#' @param ref_panel Path to reference panel WITHOUT extension (e.g., "/data/1kg_eur/EUR")
#' @param pval_threshold P-value threshold for IV selection (default: 5e-8)
#' @param clump_kb Clumping window in kb (default: 10000)
#' @param clump_r2 Clumping r² threshold (default: 0.001)
#' @param use_proxy Whether to search for proxy SNPs (default: TRUE)
#' @param proxy_r2 Minimum r² for proxy SNPs (default: 0.8)
#' @param temp_dir Temporary directory for plink output (default: tempdir())
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Harmonized data frame from TwoSampleMR::harmonise_data()
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Read GWAS data
#' exposure <- read_gwas_data("exposure.txt", ...)
#' outcome <- read_gwas_data("outcome.txt", ...)
#'
#' # Harmonize
#' harmonized <- harmonize_for_mr(
#'   exp_gwas = exposure,
#'   otc_gwas = outcome,
#'   exp_name = "BMI",
#'   otc_name = "CAD",
#'   plink_bin = "/usr/bin/plink",
#'   ref_panel = "/data/1kg_eur/EUR"
#' )
#'
#' # Extract for BMM
#' bmm_data <- extract_bmm_data(harmonized)
#' }
#' Complete GWAS Harmonization Workflow
#'
#' Performs complete harmonization matching real data evaluation:
#' 1. Read GWAS data
#' 2. Filter exposure SNPs by p-value
#' 3. Format for TwoSampleMR
#' 4. LD clumping with local reference panel
#' 5. Proxy SNP search with local plink
#' 6. Harmonization
#'
#' @param exp_gwas Exposure GWAS data frame (from read_gwas_data)
#' @param otc_gwas Outcome GWAS data frame (from read_gwas_data)
#' @param exp_name Exposure trait name
#' @param otc_name Outcome trait name
#' @param plink_bin Path to plink executable (e.g., "/usr/bin/plink")
#' @param ref_panel Path to reference panel WITHOUT extension (e.g., "/data/1kg_eur/EUR")
#' @param pval_threshold P-value threshold for IV selection (default: 5e-8)
#' @param clump_kb Clumping window in kb (default: 10000)
#' @param clump_r2 Clumping r² threshold (default: 0.001)
#' @param use_proxy Whether to search for proxy SNPs (default: TRUE)
#' @param proxy_r2 Minimum r² for proxy SNPs (default: 0.8)
#' @param temp_dir Temporary directory for plink output (default: tempdir())
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return Harmonized data frame from TwoSampleMR::harmonise_data()
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Read GWAS data
#' exposure <- read_gwas_data("exposure.txt", ...)
#' outcome <- read_gwas_data("outcome.txt", ...)
#'
#' # Harmonize
#' harmonized <- harmonize_for_mr(
#'   exp_gwas = exposure,
#'   otc_gwas = outcome,
#'   exp_name = "BMI",
#'   otc_name = "CAD",
#'   plink_bin = "/usr/bin/plink",
#'   ref_panel = "/data/1kg_eur/EUR"
#' )
#'
#' # Extract for BMM
#' bmm_data <- extract_bmm_data(harmonized)
#' }
harmonize_for_mr <- function(exp_gwas, otc_gwas,
                             exp_name = "exposure",
                             otc_name = "outcome",
                             plink_bin, ref_panel,
                             pval_threshold = 5e-8,
                             clump_kb = 10000,
                             clump_r2 = 0.001,
                             use_proxy = TRUE,
                             proxy_r2 = 0.8,
                             temp_dir = tempdir(),
                             verbose = TRUE) {
  
  # 内部警告处理函数
  precision_handler <- function(w) {
    msg <- w$message
    
    # 1. 捕获 EAF 转换/缺失警告
    if (grepl("eaf column is (not numeric|missing|invalid)", msg, ignore.case = TRUE)) {
      # 输出指定的英文消息
      message("Note: The data lacks an EAF column. The potential outlier SNP filtering step cannot be performed.")
      # 抑制原始警告
      invokeRestart("muffleWarning")
    } 
    # 2. 精确捕获 X8 列名填充警告并完全忽略
    else if (grepl("Missing column names filled in: 'X8'", msg, fixed = TRUE) ||
             grepl("X8", msg, fixed = TRUE) && grepl("filled in", msg, ignore.case = TRUE)) {
      # 完全抑制，不输出任何消息
      invokeRestart("muffleWarning")
    }
    # 3. 其他类型的警告，保持原样显示
    else {
      # 对于其他警告，保持原样显示
      # 可以在此添加其他特定的警告处理
    }
  }
  
  # 内部消息处理函数
  message_handler <- function(m) {
    msg <- m$message
    
    # 1. 捕获 "No phenotype name specified" 消息并抑制
    if (grepl("No phenotype name specified", msg, ignore.case = TRUE)) {
      # 完全抑制，不输出任何消息
      invokeRestart("muffleMessage")
    }
    # 2. 其他类型的消息，保持原样显示
    else {
      # 对于其他消息，保持原样显示
      # 可以在此添加其他特定的消息处理
    }
  }
  
  # 检查并处理名称参数
  if (missing(exp_name) || is.null(exp_name) || is.na(exp_name) || exp_name == "") {
    exp_name <- "exposure"
    if (verbose) message("Note: Using default exposure name: 'exposure'")
  }
  
  if (missing(otc_name) || is.null(otc_name) || is.na(otc_name) || otc_name == "") {
    otc_name <- "outcome"
    if (verbose) message("Note: Using default outcome name: 'outcome'")
  }
  
  if (verbose) {
    cat("\n========================================\n")
    cat(sprintf("Analyzing: %s -> %s\n", exp_name, otc_name))
    cat("========================================\n\n")
  }
  
  plink_bin <- normalizePath(plink_bin, mustWork = TRUE)
  ref_panel <- path.expand(ref_panel)
  
  # === Step 1: Validate inputs ===
  if (!file.exists(plink_bin)) {
    stop(sprintf("Plink binary not found: %s", plink_bin))
  }

  ref_panel_bed <- paste0(ref_panel, ".bed")
  if (!file.exists(ref_panel_bed)) {
    stop(sprintf("Reference panel not found: %s", ref_panel_bed))
  }

  # === Step 2: Filter exposure by p-value ===
  if (verbose) cat(sprintf("Step 2: Filtering exposure SNPs (P < %.0e)...\n", pval_threshold))

  exp_gwas_sig <- exp_gwas %>% dplyr::filter(P < pval_threshold)

  if (verbose) cat(sprintf("  %d significant SNPs\n", nrow(exp_gwas_sig)))

  if (nrow(exp_gwas_sig) < 3) {
    stop("Too few IVs (< 3)")
  }

  # === Step 3: Format for TwoSampleMR ===
  if (verbose) cat("\nStep 3: Formatting data for TwoSampleMR...\n")

  # 为 exposure 数据添加 phenotype 列
  exp_gwas_sig_with_pheno <- exp_gwas_sig %>%
    dplyr::mutate(phenotype = exp_name)

  # 对 format_data 应用警告和消息处理器
  exp_dat <- withCallingHandlers(
    {
      TwoSampleMR::format_data(
        exp_gwas_sig_with_pheno,
        type = "exposure",
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        effect_allele_col = "EA",
        other_allele_col = "OA",
        eaf_col = "EAF",
        pval_col = "P",
        chr_col = "CHR",
        pos_col = "POS",
        phenotype_col = "phenotype"
      )
    },
    warning = precision_handler,
    message = message_handler
  )

  # 确保 exposure 和 id.exposure 正确设置
  exp_dat <- exp_dat %>% dplyr::mutate(
    exposure = exp_name,
    id.exposure = exp_name
  )

  if (verbose) cat(sprintf("  Exposure formatted: %d SNPs\n", nrow(exp_dat)))

  # 为 outcome 数据添加 phenotype 列
  otc_gwas_with_pheno <- otc_gwas %>%
    dplyr::mutate(phenotype = otc_name)

  otc_dat <- withCallingHandlers(
    {
      TwoSampleMR::format_data(
        otc_gwas_with_pheno,
        type = "outcome",
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        effect_allele_col = "EA",
        other_allele_col = "OA",
        eaf_col = "EAF",
        pval_col = "P",
        chr_col = "CHR",
        pos_col = "POS",
        phenotype_col = "phenotype"
      )
    },
    warning = precision_handler,
    message = message_handler
  )

  # 确保 outcome 和 id.outcome 正确设置
  otc_dat <- otc_dat %>% dplyr::mutate(
    outcome = otc_name,
    id.outcome = otc_name
  )

  if (verbose) cat(sprintf("  Outcome formatted: %d SNPs\n", nrow(otc_dat)))

  # === Step 4: LD clumping ===
  if (verbose) cat("\nStep 4: LD clumping...\n")

  rsid_pval_exp <- exp_dat %>%
    dplyr::select(SNP, pval.exposure) %>%
    dplyr::rename(rsid = SNP, pval = pval.exposure) %>%
    dplyr::filter(!is.na(rsid), nzchar(rsid))

  rsid_pval_exp.clump <- tryCatch(
    {
      ieugwasr::ld_clump_local(
        rsid_pval_exp,
        clump_kb = clump_kb,
        clump_r2 = clump_r2,
        clump_p = 1,
        plink_bin = plink_bin,
        bfile = ref_panel
      )
    },
    error = function(e) {
      stop(sprintf("Clumping failed: %s", e$message))
    }
  )

  if (nrow(rsid_pval_exp.clump) == 0) {
    stop("No SNPs survived clumping")
  }

  exp_dat <- exp_dat %>% dplyr::filter(SNP %in% rsid_pval_exp.clump$rsid)

  if (verbose) cat(sprintf("  %d independent IVs after clumping\n", nrow(exp_dat)))

  # === Step 5: Proxy SNP search ===
  if (use_proxy) {
    if (verbose) cat("\nStep 5: Searching for proxy SNPs...\n")
    otc_dat <- get_proxy(
      exp_dat,
      otc_dat,
      r2 = proxy_r2,
      bin_plink = plink_bin,
      ref_panel = ref_panel,
      outfile = file.path(temp_dir, sprintf("%s_%s_proxy", exp_name, otc_name))
    )
  } else {
    if (verbose) cat("\nStep 5: Skipping proxy SNP search\n")
    otc_dat <- otc_dat %>% dplyr::filter(SNP %in% exp_dat$SNP)
  }

  if (verbose) cat(sprintf("  Outcome SNPs available: %d\n", nrow(otc_dat)))

  # === Step 6: Harmonization ===
  if (verbose) cat("\nStep 6: Harmonizing data...\n")

  harmonised_dat <- withCallingHandlers(
    {
      TwoSampleMR::harmonise_data(exp_dat, otc_dat)
    },
    warning = precision_handler,
    message = message_handler
  )

  if (is.null(harmonised_dat) || nrow(harmonised_dat) == 0) {
    stop("Harmonization failed")
  }

  if (verbose) {
    cat(sprintf(
      "  Harmonized: %d SNPs (mr_keep = %d)\n",
      nrow(harmonised_dat), sum(harmonised_dat$mr_keep)
    ))
    cat(sprintf("    Palindromic SNPs: %d\n", sum(harmonised_dat$palindromic, na.rm = TRUE)))
    cat(sprintf("    Ambiguous SNPs: %d\n", sum(harmonised_dat$ambiguous, na.rm = TRUE)))
    if (use_proxy) {
      cat(sprintf(
        "    Proxy SNPs used: %d\n",
        sum(!is.na(harmonised_dat$proxy.outcome) & harmonised_dat$proxy.outcome, na.rm = TRUE)
      ))
    }
    cat(sprintf("    SNPs to keep: %d\n", sum(harmonised_dat$mr_keep)))
  }

  if (sum(harmonised_dat$mr_keep) < 3) {
    stop("Too few SNPs (< 3) for MR analysis")
  }

  if (verbose) cat("\n  Analysis complete!\n\n")

  return(harmonised_dat)
}

#' Extract Data for BMM Analysis
#'
#' Extract cleaned data from harmonized output for BMM
#'
#' @param harmonised_dat Output from harmonize_for_mr()
#'
#' @return List with beta_X, beta_Y, se_X, se_Y, snps
#'
#' @export
extract_bmm_data <- function(harmonised_dat) {
  # Filter to mr_keep only
  mr_data <- harmonised_dat %>% dplyr::filter(mr_keep)

  if (nrow(mr_data) < 3) {
    stop("Fewer than 3 SNPs available for analysis")
  }

  list(
    beta_X = mr_data$beta.exposure,
    beta_Y = mr_data$beta.outcome,
    se_X = mr_data$se.exposure,
    se_Y = mr_data$se.outcome,
    snps = mr_data$SNP,
    n_snps = nrow(mr_data)
  )
}
