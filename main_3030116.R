# ════════════════════════════════════════════════════════════════════════════════
# PRS-MIGRAINE ASSOCIATION ANALYSIS PIPELINE ----
# ════════════════════════════════════════════════════════════════════════════════
# DESCRIPTION: Comprehensive analysis pipeline for testing polygenic risk score (PRS) 
#              associations with migraine using UK Biobank data
# AUTHORS: Amparo Gimenez Rios
# DATE: August 2025
#
# USAGE: Run source("D:/Uni/MASTER/Classes/MSc Project/Data/main_3030116.R") in console
# ════════════════════════════════════════════════════════════════════════════════

# ────────────────────────────────────────────────────────────────────────────────
# 0) PROMPTS — ask everything up front ----
# ────────────────────────────────────────────────────────────────────────────────
# PURPOSE: Collect all analysis parameters through interactive prompts
# WHY: Ensures reproducible analysis while allowing flexible parameter specification

message("To check function documentation see functions_documentation.R script.")

#region INTERACTIVE SETUP

# 0. Run identification ----
run_label <- readline("0) Enter a short label for this run (e.g. base, demog, full): ")
# in my case: base (PGC8), demog (sex, age, townsend, chip and PGC8) and full (all)

# Analysis scope configuration ----
# 1. ask if preliminary exploration should run
prelim_input <- readline("1) Run preliminary EDA overview (demographics & tests)? (y/n, default y): ")
do_prelim <- !nzchar(prelim_input) || tolower(trimws(prelim_input)) %in% c("y","yes")

# 2. ask if full EDA should run
do_eda_input <- readline("2) Do you want to run EDA plots now? (y/n, default n): ")
do_eda <- tolower(trimws(do_eda_input)) %in% c("y","yes")

# 2.1. EDA variable specification ----
# if preliminary or full EDA are specified, ask for relevant variables
if (do_prelim || do_eda) {
  message("EDA variables need to be specified.")
  age_var <- readline("2.1) Age variable column name (e.g. age_0_0): ")
  chip_var <- readline("2.2) Chip variable column name (e.g. chip): ")
  fem_level <- readline("2.3) Specify female level in sex variable): ")

} else {
  message("Skipping preliminary EDA and full EDA plots.")
  age_var  <- ""
  chip_var <- ""
  fem_level <- ""
}
# age_0_0
# chip
# 2

# 3. Core model variables ----
# sex variable is needed for the models too, so ask for it here outside the if statement
sex_var <- readline("3) Sex variable column name (e.g. sex_plink): ")
# sex_plink

# Validate that sex variable is provided since it's required for modeling
if (!nzchar(sex_var)) {
  stop("❌ Sex variable is required for PRS association modeling. Please provide a valid column name.")
}

# 4. File paths and data specification ----
message("Input directory must contain data and R scripts (functions and main)")
input_dir  <- readline("4) Full path to your DATA directory: ")
if (!nzchar(input_dir) ||
    !grepl("^(?:[A-Za-z]:[\\\\/]|/|~)", input_dir, perl = TRUE)) {
  stop("✗ Invalid path supplied for data directory: ", input_dir)
}
if (!dir.exists(input_dir)) {
  stop("✗ Data directory not found: ", input_dir)
}
# D:/Uni/MASTER/Classes/MSc Project/Data

# 5. Output path
message("Output directory will be created if it does not exist.")
output_dir <- readline("5) Full path to your RESULTS directory: ")
if (!nzchar(output_dir) ||
    !grepl("^(?:[A-Za-z]:[\\\\/]|/|~)", output_dir, perl = TRUE)) {
  stop("✗ Invalid path supplied for results directory: ", output_dir)
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
# D:/Uni/MASTER/Classes/MSc Project/Data/Results

# 6. UK BB dataframe name
ukb_fn     <- readline("6) UKB .dta filename (in your data dir): ")
# UKB_71392_migraine_project.dta
# ukbb_subset.txt

# 7. Variable column specification ----
id_col        <- readline("7) ID column name (e.g. n_eid): ")
# n_eid

# 8. Dependent variable
status_col    <- readline("8) Dependent variable status column name (levels 0/1): ")
# n_120016_0_0

# 9. Independent variables
prs_cols      <- strsplit(
  readline("9) Independent variables PRS column names (comma‐sep, e.g. prs1,prs2): "),
  ",\\s*"
)[[1]]
# prs_asd_std, prs_adhd_std

# 10. Covariables
covariates    <- strsplit(
  readline("10) Covariate column names (comma‐sep, or blank): "),
  ",\\s*"
)[[1]]
# base: PGC1, PGC2, PGC3, PGC4, PGC5, PGC6, PGC7, PGC8
# demog: sex_plink, age_0_0, townsend, PGC1, PGC2, PGC3, PGC4, PGC5, PGC6, PGC7, PGC8, chip
# full: sex_plink, age_0_0, townsend, PGC1, PGC2, PGC3, PGC4, PGC5, PGC6, PGC7, PGC8, chip, n_21001_0_0, curr_alcohol_0, cursmoke_0
# (includes antropometric, lifestyle, and other health variables)

# Statistical analysis parameters ----
message("The next prompts can be skipped to use default values.")
# 11. VIF threshold
vif_input <- readline("11) VIF threshold (default 5 - enter) for multicollinearity: ")
vif_thr   <- if (nzchar(vif_input)) as.numeric(vif_input) else 5
# if default is desired, press enter when prompted without writing anything

# 12. Hosmer-Lemeshow groups threshold
hl_input  <- readline("11) Hosmer–Lemeshow groups (default 10 - enter) to compute χ² goodness-of-fit: ")
hl_groups <- if (nzchar(hl_input)) as.integer(hl_input) else 10
# if default is desired, press enter when prompted without writing anything

# 13. Cooks assumption threshold
cooks_value <- readline("12) Set cutoff value for Cooks assumption (if not specified, 1 for large datasets): ")
cooks_cutoff <- if (nzchar(cooks_value)) as.integer(cooks_value) else NULL

if (nzchar(cooks_value)) {
  cooks_cutoff <- as.numeric(cooks_value)
  if (is.na(cooks_cutoff) || cooks_cutoff <= 0) {
    warning("Invalid cut‑off supplied; using default 1 instead.")
    cooks_cutoff <- NULL
  }
} else {
  cooks_cutoff <- NULL  # triggers default in cooks_check()
}

# Path validation and setup ----
# Sanity checks / build paths
if (!dir.exists(input_dir))  stop("Data directory not found: ", input_dir)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ukb_path     <- file.path(input_dir, ukb_fn)
#endregion


# ────────────────────────────────────────────────────────────────────────────────
# 1. Load functions and packages ----
# ────────────────────────────────────────────────────────────────────────────────
# PURPOSE: Import all required R functions and packages for the PRS analysis pipeline
# WHY: We need custom functions for data processing, modeling, and assumption checking
#      as well as standard R packages for statistical analysis and visualization

#region SETUP FUNCTIONS AND PACKAGES
functions_path <- file.path(input_dir, "functions_3030116.R")
# this saves the complete path to functions_project.R
if (! file.exists(functions_path)) {
  # if the file does not exist, the script stops and gives an error message
  stop("Could not find functions_3030116.R in: ", input_dir)
}

source(functions_path)  # Load all custom functions from external file

# install/load all packages needed
load_packages()  # Custom function that loads required packages (ggplot2, pROC, etc.)
set_plot_defaults(base_size = 16)  # Set a readable global ggplot theme
# packages needed for MVMR are not included and will be loaded later if pipeline
# proceeds to MVMR
#endregion


# ────────────────────────────────────────────────────────────────────────────────
# 2. Load and clean UK Biobank data ----
# ────────────────────────────────────────────────────────────────────────────────
# PURPOSE: Load raw UK Biobank data and perform initial cleaning steps
# WHY: Raw data needs cleaning (missing values, outliers) and standardization
#      before statistical modeling can be performed safely and reliably

#region DATA LOADING AND CLEANING
ukbb       <- load_ukbb_data(ukb_path)  # Load raw UK Biobank dataset
ukbb       <- clean_migraine_status(ukbb, migraine_col = status_col)  # Clean outcome variable

# Validate that all requested columns exist in the dataset
all_requested_vars <- unique(c(id_col, status_col, prs_cols, sex_var,
                              if(nzchar(age_var)) age_var else character(0),
                              if(nzchar(chip_var)) chip_var else character(0),
                              covariates))
missing_vars <- setdiff(all_requested_vars, names(ukbb))

if (length(missing_vars) > 0) {
  cat("\n❌ ERROR: The following variables are not found in your dataset:\n")
  cat(paste("-", missing_vars, collapse = "\n"), "\n")
  cat("\nAvailable columns in your dataset:\n")
  cat(paste("-", sort(names(ukbb)), collapse = "\n"), "\n")
  stop("Please check your variable names and try again.")
}

# to select columns, it is needed to take into account that if EDA is in place and covariates are used,
# the variables may overlap, so only unique values should be kept
add_if <- function(x) {
  if (length(x) == 0) return(character(0))
  x[nzchar(x)]
}
# add_if will check if the variable is not empty, and if it is not, it will return it
eda_vars <- unique(c(sex_var, add_if(age_var), add_if(chip_var), add_if(covariates)))

ukbb_clean <- select_columns(
  ukbb,
  cols_to_keep = unique(c(id_col, status_col, prs_cols, eda_vars))  # Keep only necessary columns
)


ukbb_clean <- remove_missing_prs(
  ukbb_clean,
  cols_na = prs_cols  # Remove individuals with missing PRS values
)

# Normalize binary covariates (e.g., smoking/alcohol) to 0/1 while preserving NAs
# This ensures glm treats them as proper binary indicators even if source uses 0/2 coding
to_bin01 <- function(x) {
  if (is.numeric(x)) return(ifelse(is.na(x), NA, as.integer(x != 0)))
  x
}
if ("cursmoke_0" %in% names(ukbb_clean)) {
  ukbb_clean$cursmoke_0 <- to_bin01(ukbb_clean$cursmoke_0)
}
if ("curr_alcohol_0" %in% names(ukbb_clean)) {
  ukbb_clean$curr_alcohol_0 <- to_bin01(ukbb_clean$curr_alcohol_0)
}

ukbb_no <- remove_outliers_z(
  df           = ukbb_clean,
  numeric_vars = prs_cols  # Create version without extreme PRS outliers (>5 SD)
)
# Modeling-only median imputation for BMI (keeps EDA raw)
impute_median_var <- function(df, var) {
  if (var %in% names(df)) {
    med <- suppressWarnings(median(df[[var]], na.rm = TRUE))
    if (is.finite(med)) {
      n_before <- sum(is.na(df[[var]]))
      if (n_before > 0) {
        df[[var]][is.na(df[[var]])] <- med
        message(sprintf("Imputed %d NA(s) in %s with median %.3f (modeling only)", n_before, var, med))
      }
    }
  }
  df
}

# Create modeling copies with BMI imputed (if present)
ukbb_clean_model <- impute_median_var(ukbb_clean, "n_21001_0_0")
ukbb_no_model    <- impute_median_var(ukbb_no,    "n_21001_0_0")
#endregion


# ────────────────────────────────────────────────────────────────────────────────
# 3. Optional preliminary exploration ----
# ────────────────────────────────────────────────────────────────────────────────
# PURPOSE: Generate descriptive statistics and basic exploratory data analysis
# WHY: Understanding data characteristics before modeling helps identify issues
#      and provides context for interpreting results

#region PRELIMINARY EXPLORATION
#ensure female level is not empty
if (do_prelim && (!nzchar(fem_level))) {
  fem_level <- readline(sprintf("What value in '%s' indicates female? ", sex_var))  # Interactive prompt for sex coding
}

# proceed to preliminary EDA if requested
if (do_prelim) {
  
  prelim <- prelimin_eda(
    df         = ukbb_no,
    status_var = status_col,
    sex_var    = sex_var,
    age_var    = age_var,
    chip_var   = chip_var,
    female_code = fem_level,
    output_dir = output_dir,
    outfile      = "prelimin_eda_summary.txt", # Generate comprehensive EDA report
  dataset_label = "without_outliers",
  n_with_outliers = nrow(ukbb_clean)
  )

    prelim <- prelimin_eda(
    df         = ukbb_clean,
    status_var = status_col,
    sex_var    = sex_var,
    age_var    = age_var,
    chip_var   = chip_var,
    female_code = fem_level,
    output_dir = output_dir,
    outfile      = "prelimin_eda_summary.txt", # Generate comprehensive EDA report
  dataset_label = "with_outliers"
  )
  
  message("Preliminary EDA written to: ", file.path(output_dir, "prelimin_eda_summary.txt"))
}

# EDA plots
if (do_eda) {
  run_all_eda(
    ukbb_clean = ukbb_clean,
    prs_cols   = prs_cols,
    facet_vars = c(sex_var, chip_var),
    status_col = status_col,
    covariates = covariates,
    output_dir = output_dir  # Generate violin plots, decile plots, correlation heatmaps, etc.
  )
}
#endregion


# ────────────────────────────────────────────────────────────────────────────────
# 4. PRS association modeling and assumption checking ----
# ────────────────────────────────────────────────────────────────────────────────
# PURPOSE: Fit multiple logistic regression models to test PRS-migraine associations
# WHY: We test different interaction patterns to understand how PRS effects vary
#      by context, and we need to verify that statistical assumptions are met

#region PRS ASSOCIATION MODELING
message("Running association analyses and assumption checks…")

data_list <- list(
  with_outliers    = ukbb_clean_model,    # Use BMI-imputed for modeling
  without_outliers = ukbb_no_model        # Imputed version of outlier-removed dataset
)

# Remove sex from covariates for stratified runs (no interactions; sex not included)
covariates_nosex <- setdiff(covariates, sex_var)

# directory for assumption plots
ass_dir <- file.path(output_dir, "Assumptions")
if (!dir.exists(ass_dir)) dir.create(ass_dir, recursive = TRUE)

# initialise empty objects to store results
# Initialize output objects for models (with and without outliers)
out_no_int_with      <- NULL  # No interaction model (with outliers)
out_prsprs_with      <- NULL  # PRS×PRS interaction model (with outliers)
out_sex_with         <- NULL  # PRS×Sex interaction model (with outliers)
out_both_with        <- NULL  # Both interactions model (with outliers)

out_no_int_without   <- NULL  # No interaction model (without outliers)
out_prsprs_without   <- NULL  # PRS×PRS interaction model (without outliers)
out_sex_without      <- NULL  # PRS×Sex interaction model (without outliers)
out_both_without      <- NULL  # Both interactions model (without outliers)

# Initialize output objects for assumption checks (with and without outliers)
chk_no_int_with      <- NULL  # Assumption checks for no interaction (with outliers)
chk_prsprs_with      <- NULL  # Assumption checks for PRS×PRS interaction (with outliers)
chk_sex_with         <- NULL  # Assumption checks for PRS×Sex interaction (with outliers)
chk_both_with        <- NULL  # Assumption checks for both interactions (with outliers)

chk_no_int_without   <- NULL  # Assumption checks for no interaction (without outliers)
chk_prsprs_without   <- NULL  # Assumption checks for PRS×PRS interaction (without outliers)
chk_sex_without      <- NULL  # Assumption checks for PRS×Sex interaction (without outliers)
chk_both_without      <- NULL  # Assumption checks for both interactions (without outliers)

# Initialize containers for sex-stratified results and their assumption checks
stratified_with            <- NULL # Sex-stratified results (with outliers)
stratified_without         <- NULL # Sex-stratified results (without outliers)
stratified_checks_with     <- NULL # Assumption checks for sex-stratified models (with outliers)
stratified_checks_without  <- NULL # Assumption checks for sex-stratified models (without outliers)


# loop over each version
for (name in names(data_list)) {
  df <- data_list[[name]]
  
  ## a) fit models ---------------------------------------------------------------
  # We fit 4 different models to test different interaction hypotheses:
  
  # 1) No interaction at all - tests main effects of each PRS independently
  out_no_int <- fit_prs_model(
    df, outcome_var=status_col, prs_vars=prs_cols,
    covariates=covariates,
    interaction=FALSE,           # No PRS×PRS interaction
    sex_var=sex_var,           # you can supply sex_var even if sex_interaction=FALSE
    sex_interaction=FALSE        # No PRS×Sex interaction
  )
  chk_no_int <- check_all_assumptions(
    fit=out_no_int$model, numeric_vars=prs_cols,
    vif_thresh=vif_thr, vif_type="terms",      # Check multicollinearity
    hl_groups=hl_groups, cook_cutoff=cooks_cutoff  # Check goodness of fit and outliers
  )
  
  # 2) PRS × PRS only - tests whether ADHD and ASD PRS effects interact
  out_prsprs <- fit_prs_model(
    df, outcome_var=status_col, prs_vars=prs_cols,
    covariates=covariates,
    interaction=TRUE,            # Include PRS×PRS interaction
    sex_var=sex_var,
    sex_interaction=FALSE        # No PRS×Sex interaction
  )
  chk_prsprs <- check_all_assumptions(
    fit=out_prsprs$model, numeric_vars=prs_cols,
    vif_thresh=vif_thr, vif_type="predictor",  # Use stricter VIF check for interaction model
    hl_groups=hl_groups, cook_cutoff=cooks_cutoff
  )
  
  # 3) PRS × Sex only - tests whether PRS effects differ by biological sex
  out_sex <- fit_prs_model(
    df, outcome_var=status_col, prs_vars=prs_cols,
    covariates=covariates,
    interaction=FALSE,           # No PRS×PRS interaction
    sex_var=sex_var,
    sex_interaction=TRUE         # Include PRS×Sex interaction
  )
  chk_sex <- check_all_assumptions(
    fit=out_sex$model, numeric_vars=prs_cols,
    vif_thresh=vif_thr, vif_type="terms",
    hl_groups=hl_groups, cook_cutoff=cooks_cutoff
  )
  
  # 4) Both interactions together - most complex model with all interactions
  out_both <- fit_prs_model(
    df, outcome_var=status_col, prs_vars=prs_cols,
    covariates=covariates,
    interaction=TRUE,            # Include PRS×PRS interaction
    sex_var=sex_var,
    sex_interaction=TRUE         # Include PRS×Sex interaction
  )
  chk_both <- check_all_assumptions(
    fit=out_both$model, numeric_vars=prs_cols,
    vif_thresh=vif_thr, vif_type="predictor",  # Use stricter VIF check for full interaction model
    hl_groups=hl_groups, cook_cutoff=cooks_cutoff
  )

  # 5) Stratified by sex (no interactions; sex removed from covariates)
  strat_res <- run_sex_stratified_no_interactions(
    df = df,
    outcome_var = status_col,
    prs_vars = prs_cols,
    covariates = covariates_nosex,
    sex_var = sex_var
  )
  
  # Assumption checks and plots for stratified (no-interaction) models per sex level
  strat_checks <- list()
  if (!is.null(strat_res) && length(strat_res) > 0) {
    for (lv in names(strat_res)) {
      sres <- strat_res[[lv]]
      if (is.null(sres) || is.null(sres$model)) next
      chk_s <- check_all_assumptions(
        fit = sres$model,
        numeric_vars = prs_cols,
        vif_thresh = vif_thr,
        vif_type = "terms",
        hl_groups = hl_groups,
        cook_cutoff = cooks_cutoff
      )
      strat_checks[[lv]] <- chk_s
      # Save plots for each sex stratum
      safe_lv <- gsub("[^A-Za-z0-9_]+", "_", as.character(lv))
      suf <- paste0("_noInt_strat_", name, "_sex_", safe_lv)
      save_crs(chk_s, suf, ass_dir)
      ggsave(file.path(ass_dir, paste0("roc", suf, ".png")),
             chk_s$roc_auc, width = 10, height = 8, dpi = 300)
      ggsave(file.path(ass_dir, paste0("independence", suf, ".png")),
             chk_s$independence, width = 10, height = 8, dpi = 300)
    }
  }
  
  ## b) store into with/without buckets -----------------------------------------
  if (name == "with_outliers") {
    out_no_int_with    <- out_no_int
    out_prsprs_with    <- out_prsprs
    out_sex_with       <- out_sex
    out_both_with      <- out_both
    
    chk_no_int_with    <- chk_no_int
    chk_prsprs_with    <- chk_prsprs
    chk_sex_with       <- chk_sex
    chk_both_with      <- chk_both
  stratified_with    <- strat_res
  stratified_checks_with <- strat_checks
  } else {
    out_no_int_without <- out_no_int
    out_prsprs_without <- out_prsprs
    out_sex_without    <- out_sex
    out_both_without   <- out_both
    
    chk_no_int_without <- chk_no_int
    chk_prsprs_without <- chk_prsprs
    chk_sex_without    <- chk_sex
    chk_both_without   <- chk_both
  stratified_without <- strat_res
  stratified_checks_without <- strat_checks
  }
  
  ## c) Generate statistical assumption diagnostic plots ----------------------
  # Create comprehensive diagnostic visualizations for each model to validate
  # the underlying statistical assumptions of logistic regression
  
  models <- list(
    noInt    = list(chk=chk_no_int,    suf="_noInt",    name="No Interactions"),
    prsprs   = list(chk=chk_prsprs,    suf="_prsprs",   name="PRS × PRS Only"),
    sex      = list(chk=chk_sex,       suf="_sex",      name="PRS × Sex Only"),
    both     = list(chk=chk_both,      suf="_both",     name="Both Interactions")
  )
  
  for (m in names(models)) {
    suf <- models[[m]]$suf
    chk <- models[[m]]$chk
    model_name <- models[[m]]$name
    
    # Component residual plots - check linearity assumption
    # Shows whether PRS have linear relationships with log-odds
    save_crs(chk, paste0(suf,"_",name), ass_dir)
    
    # ROC curves - visualize model discrimination ability
    # AUC values indicate how well models separate cases from controls
  save_plot_hq(file.path(ass_dir,paste0("roc",suf,"_",name,".png")),
         chk$roc_auc, width=10, height=8, dpi=300)
    
    # Independence plots - check residual patterns  
    # Random scatter suggests independence assumption is met
  save_plot_hq(file.path(ass_dir,paste0("independence",suf,"_",name,".png")),
         chk$independence, width=10, height=8, dpi=300)
  }
}

message("Saved assumption plots for all interaction types (noInt/prsprs/sex/both) and both data versions.")

# ────────────────────────────────────────────────────────────────────────────────
# 5. Model results compilation and export ----  
# ────────────────────────────────────────────────────────────────────────────────
# PURPOSE: Compile comprehensive results from all models into exportable format
# WHY: Create standardized output tables for statistical reporting and comparison
#      across different interaction models and outlier handling approaches

#region MODEL RESULTS EXPORT
# Organize all fitted models for results compilation
# Each model represents a different hypothesis about PRS interactions
model_list <- list(
  noInt_with    = out_no_int_with$model,     # Main effects only (with outliers)
  prsprs_with   = out_prsprs_with$model,     # PRS×PRS interaction (with outliers)
  sex_with      = out_sex_with$model,        # PRS×Sex interaction (with outliers)  
  both_with     = out_both_with$model,       # Both interactions (with outliers)
  noInt_without = out_no_int_without$model,  # Main effects only (without outliers)
  prsprs_without= out_prsprs_without$model,  # PRS×PRS interaction (without outliers)
  sex_without   = out_sex_without$model,     # PRS×Sex interaction (without outliers)
  both_without  = out_both_without$model     # Both interactions (without outliers)
)

# Organize all assumption check results for comprehensive reporting
# Each check includes VIF, Hosmer-Lemeshow test, Cook's distance, and plots
check_list <- list(
  noInt_with    = chk_no_int_with,     # Assumption checks for main effects (with outliers)
  prsprs_with   = chk_prsprs_with,     # Assumption checks for PRS×PRS (with outliers)
  sex_with      = chk_sex_with,        # Assumption checks for PRS×Sex (with outliers)
  both_with     = chk_both_with,       # Assumption checks for both interactions (with outliers)
  noInt_without = chk_no_int_without,  # Assumption checks for main effects (without outliers)
  prsprs_without= chk_prsprs_without,  # Assumption checks for PRS×PRS (without outliers)
  sex_without   = chk_sex_without,     # Assumption checks for PRS×Sex (without outliers)
  both_without  = chk_both_without     # Assumption checks for both interactions (without outliers)
)

# Export consolidated results table with all key statistics:
# - Model coefficients (OR, 95% CI, p-values) for each PRS term
# - Model fit statistics (AUC, Nagelkerke R², AIC)
# - Assumption test results (Hosmer-Lemeshow p-value, max VIF)
# - Outlier diagnostics (Cook's distance flags)
export_consolidated_results(
  model_list = model_list,
  check_list = check_list,
  output_dir = output_dir,
  run_label  = run_label
)


message("Saved combined model term details to CSV.")

# Export sex-stratified (no-interaction) summaries
make_strat_df <- function(strat_list, dataset_label) {
  if (is.null(strat_list) || length(strat_list) == 0) return(NULL)
  lvls <- names(strat_list)
  parts <- lapply(lvls, function(lv) {
    res <- strat_list[[lv]]
    if (is.null(res)) return(NULL)
    df <- res$results
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$sex_level <- lv
    df
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  out <- do.call(rbind, parts)
  out$dataset <- dataset_label
  out
}

strat_with_df    <- make_strat_df(stratified_with,    "with_outliers")
strat_without_df <- make_strat_df(stratified_without, "without_outliers")

if (!is.null(strat_with_df)) {
  write.csv(
    strat_with_df,
    file = file.path(output_dir, paste0("results_stratified_sex_", run_label, "_with.csv")),
    row.names = FALSE
  )
}
if (!is.null(strat_without_df)) {
  write.csv(
    strat_without_df,
    file = file.path(output_dir, paste0("results_stratified_sex_", run_label, "_without.csv")),
    row.names = FALSE
  )
}

# Export diagnostics summary for sex-stratified models
make_strat_diag_df <- function(strat_list, check_list, dataset_label) {
  if (is.null(strat_list) || is.null(check_list)) return(NULL)
  lvls <- intersect(names(strat_list), names(check_list))
  if (!length(lvls)) return(NULL)
  parts <- lapply(lvls, function(lv) {
    sres <- strat_list[[lv]]
    chk  <- check_list[[lv]]
    if (is.null(sres) || is.null(chk) || is.null(sres$model)) return(NULL)
    mod <- sres$model
    y     <- stats::model.response(mod$model)
    preds <- stats::fitted(mod)
    auc_v <- tryCatch(as.numeric(pROC::auc(y, preds)), error = function(e) NA_real_)
    data.frame(
      dataset = dataset_label,
      sex_level = lv,
      AIC = AIC(mod),
  Nagelkerke_R2 = tryCatch(as.numeric(pscl::pR2(mod)["r2CU"]), error = function(e) NA_real_),
      AUC = auc_v,
      VIF_pass = chk$multicol$pass,
      VIF_max = suppressWarnings(max(chk$multicol$vif_values, na.rm = TRUE)),
      Cooks_pass = chk$influencers$pass,
      Cooks_flagged = chk$influencers$n_flagged,
      Hosmer_Lemeshow_p = chk$hosmer$p.value,
      HL_pass = chk$hosmer$pass,
      stringsAsFactors = FALSE
    )
  })
  parts <- Filter(Negate(is.null), parts)
  if (!length(parts)) return(NULL)
  do.call(rbind, parts)
}

diag_with_df    <- make_strat_diag_df(stratified_with,    stratified_checks_with,    "with_outliers")
diag_without_df <- make_strat_diag_df(stratified_without, stratified_checks_without, "without_outliers")

diag_combined <- rbind(
  if (!is.null(diag_with_df)) diag_with_df else NULL,
  if (!is.null(diag_without_df)) diag_without_df else NULL
)
if (!is.null(diag_combined) && nrow(diag_combined) > 0) {
  write.csv(
    diag_combined,
    file = file.path(output_dir, paste0("results_stratified_sex_", run_label, "_diagnostics.csv")),
    row.names = FALSE
  )
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 6. Cross-validation performance assessment ----  
# ────────────────────────────────────────────────────────────────────────────────
#region CROSS-VALIDATION
# Perform 5-fold cross-validation to assess model generalizability and 
# obtain unbiased estimate of predictive performance using main effects model

# 1) Build formula for main effects model (most parsimonious)
formula_str <- paste(
  status_col, 
  "~", 
  paste(c(prs_cols, covariates), collapse = " + ")
)

# 2) Run 5-fold cross-validation with reproducible seed
set.seed(123)
res_cvAUC <- cv_auc(
  df          = ukbb_clean_model,      # Use full dataset including outliers (BMI-imputed)
  formula_str = formula_str,     # Main effects model formula
  k           = 5,               # Number of folds
  seed        = 123              # Reproducibility seed
)

# 3) Calculate summary statistics and save results  
mean_auc <- mean(res_cvAUC$aucs, na.rm = TRUE)    # Average AUC across folds
se_auc   <- sd(res_cvAUC$aucs, na.rm = TRUE) / sqrt(length(res_cvAUC$aucs))  # Standard error

# Create results dataframe for export
df_cv <- data.frame(
  method    = "cvAUC",
  estimate  = mean_auc,
  std_error = se_auc
)

# Export cross-validation results to CSV
write.csv(
  df_cv,
  file = file.path(output_dir, "cvAUC_results.csv"),
  row.names = FALSE
)
message("cvAUC results saved to cvAUC_results.csv")

# 4) Create visualization of cross-validation performance
# Generate bar plot showing AUC estimate with 95% confidence intervals
p_cv <- ggplot(df_cv, aes(x = method, y = estimate, fill = method)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(
    ymin = estimate - 1.96 * std_error,    # Lower 95% CI
    ymax = estimate + 1.96 * std_error     # Upper 95% CI
  ), width = 0.15) +
  labs(x = NULL, y = "AUC") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 16),
    axis.text.y  = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16)
  ) +
  coord_cartesian(ylim = c(0, 1))

# Save cross-validation performance plot
save_plot_hq(
  file.path(output_dir, "cvAUC_barplot.png"),
  p_cv, width = 6, height = 4, dpi = 300
)
message("Bar plot saved to cvAUC_barplot.png")
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 7. Pipeline completion and model archiving ----
# ────────────────────────────────────────────────────────────────────────────────
# PURPOSE: Save all fitted models and notify user of successful completion
# WHY: Preserve model objects for future analysis and provide clear completion status

#region PIPELINE COMPLETION
# Print success message to inform user of completed analysis
message("🎉 Pipeline ran successfully!")

# Archive all fitted model objects in a single RData file for this analysis run
# This preserves the complete set of models (with/without outliers, all interaction types)
#endregion
