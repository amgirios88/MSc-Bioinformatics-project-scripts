# functions.R

message("functions_3030116.R script successfully loaded.") 

# ────────────────────────────────────────────────────────────────────────────────
# Script structure
# ────────────────────────────────────────────────────────────────────────────────
# Helper functions for PRS ⇒ migraine association workflow
# (data I/O, cleaning, EDA, modeling, diagnostics, exports, CV)
#
# Sections and function order:
# 1) Setup and I/O (Input/Output)
#    - load_packages()
#    - load_ukbb_data()
#    - clean_migraine_status()
#    - select_columns()
#    - remove_missing_prs()
#    - remove_outliers_z()
#
# 2) Small helpers (internal)
#    - .ensure_factor01()
#    - .infer_case_level()
#    - .as_binary01()
#    - .predict_case_prob()
        # Note: the dot convention was suggested by Copilot. I am unsure if this is used since I couldn't find much information about it, but I decided to
        # apply it anyway because I think it is a good way to differentiate internal helper functions from regular functions.
#
# 3) Plot themes
#    - plot_theme_readable()
#
# 4) EDA plots and utilities
#    - plot_sex_distribution()
#    - plot_chip_distribution()
#    - choose_2sample_test()
#    - prelimin_eda()
#    - plot_prs_decile()
#    - plot_uni_or()
#    - eda_save_prs_decile()
#    - eda_save_uni_or()
#    - run_all_eda()
#
# 5) Modeling
#    - fit_prs_model()
#    - run_sex_stratified_no_interactions()
#
# 6) Assumption checks and diagnostics
#    - plot_cr_gg()
#    - vif_check()
#    - cooks_check()
#    - hoslem_check()
#    - plot_resid_index()
#    - check_all_assumptions()
#
# 7) Utilities for saving/exports
#    - set_plot_defaults()
#    - save_plot_hq()
#    - save_crs()
#    - export_consolidated_results()
#
# 8) Cross‑validation and ROC
#    - diagnose_cv_data()
#    - cv_auc()
#    - plot_roc_gg()
# ────────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────────
# 1. Load all required packages ----

# this function checks if required packages are installed, installs them if not,
# and loads them into the R session, with a tryCatch block to handle errors

#region LOAD PACKAGES
load_packages <- function() {
  message("Checking required packages …")
  
  # list of required packages
  required <- c(
    # data I/O and manipulation
  "haven", "dplyr", "ggplot2", "broom", "purrr", "tibble",
    # modeling and diagnostics
    "car", "ResourceSelection", "pROC", "caret", "rms", "pscl"
  )
  
  # check if each package is installed, install if not
  for(pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      # requireNamespace() checks if the package is installed without loading it
      message("Installing missing package: ", pkg)
      
      tryCatch( # tryCatch(expr, error = handler) 
        # expr is the expression to evaluate, and error is the handler function that 
        # will be called if an error occurs
        install.packages(pkg), # tries to install the package
        error = function(e) stop("Failed to install ", pkg, ": ", e$message)
      ) 
      # R creates an error condition object and hands it to error = handler
      # handler must be a function, so we define it as an anonymous function(e)
      # that takes the error object e and stops execution with a message
      #e$message accesses the error message from the error object
      
    }
    # now load pre-existing or installed packages
    tryCatch( # tryCatch(expr, error = handler) 
      library(pkg, character.only = TRUE), # expr
      error = function(e) stop("Failed to load ", pkg, ": ", e$message)
      # handler
    )
  }
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 2. Load UKB .dta ----

# this function reads a UK Biobank .dta, .csv, or .txt file using the haven package
# this is made for the UK Biobank data that was provided by supervisor, 
# which is a Stata file, with a tryCatch block to handle errors

#region LOAD UKB DATA
load_ukbb_data <- function(filepath) {
  if (!file.exists(filepath)) 
    stop("File not found: ", filepath)
  
  ext <- tolower(tools::file_ext(filepath))
  
  df <- switch(ext,
               # switch() tests an expression against elements of a list. If the value 
               # evaluated from the expression matches an item from the list, the 
               # corresponding value is returned
               # https://www.datamentor.io/r-programming/switch-function
               
               # if the file extension is .dta, read it using haven::read_dta
               dta = tryCatch(
                 haven::read_dta(filepath),
                 error = function(e) stop("Error reading .dta: ", e$message)
               ),
               
               # if the file extension is .txt or .csv read it using read.csv(..., sep="")
               txt = tryCatch(
                 read.csv(filepath, sep = "", header = TRUE, stringsAsFactors = FALSE),
                 error = function(e) stop("Error reading .txt: ", e$message)
               ),
               csv = tryCatch(
                 read.csv(filepath, header = TRUE, stringsAsFactors = FALSE),
                 error = function(e) stop("Error reading .csv: ", e$message)
               ),
               stop("Unsupported file extension: .", ext)
  )
  
  df
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 3. Clean migraine status ----

# this function filters the dataframe to keep only rows where the migraine status
# column (default 'n_120016_0_0') is either 0 or 1, and removes rows with NA values, 
# controls for errors

#region CLEAN MIGRAINE STATUS
clean_migraine_status <- function(df, migraine_col) {
  if (!is.data.frame(df)) stop("`df` must be a data.frame")
  # if df is not a data frame, stop execution and return an error message
  
  if (!migraine_col %in% names(df)) stop("Column not found: ", migraine_col)
  # if the specified migraine column is not in the data frame, stop execution and 
  # return an error message
  
  out <- df %>%
    filter(!is.na(.data[[migraine_col]]), .data[[migraine_col]] %in% c(0,1))
  # filter the data frame to keep only rows where the migraine column is not NA
  
  if (nrow(out) == 0) warning("No rows with migraine status 0 or 1")
  # if the filtered data frame is empty, issue a warning
  
  out # return the filtered data frame with only rows where migraine status is 0 or 1
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 4. Select relevant columns ----

# this function selects only the columns that are needed for the PRS‐migraine analysis
# allows for user to specify which columns to keep, controls for errors

#region SELECT COLUMNS
select_columns <- function(df, cols_to_keep) {
  # cols_to_keep is a character vector of column names to retain
  
  if (!is.data.frame(df)) stop("`df` must be a data.frame")
  # if df is not a data frame, stop execution and return an error message
  
  if (length(cols_to_keep) == 0) 
    stop("You must supply at least one column name in `cols_to_keep`")
  # if cols_to_keep is empty, stop execution and return an error message
  
  missing <- setdiff(cols_to_keep, names(df))
  # find columns that are in cols_to_keep but not in df
  
  if (length(missing)) stop("Columns not found: ", paste(missing, collapse = ", "))
  # if there are any missing columns, stop execution and return an error message
  
  out <- df[, cols_to_keep, drop = FALSE] # select only the specified columns
  
  out # return the data frame with only the specified columns
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 4. Drop missing PRS ----

# this function removes rows from the data frame where any of the specified PRS columns
# contain NA values. It is designed to work with PRS data;
# it allows the user to specify which columns to check for NA values
# and controls for missing columns

#region REMOVE MISSING PRS
remove_missing_prs <- function(df, cols_na) {
  if (!is.data.frame(df)) stop("`df` must be a data.frame")
  # if df is not a data frame, stop execution and return an error message
  
  if (length(cols_na) == 0) 
    stop("You must supply at least one column name in `cols_na`")
  # if cols_na is empty, stop execution and return an error message
  
  missing <- setdiff(cols_na, names(df))
  # find columns that are in cols_na but not in df
  
  if (length(missing)) stop("Columns not found: ", paste(missing, collapse = ", "))
  # if there are any missing columns, stop execution and return an error message
  
  out <- df %>% filter(if_all(all_of(cols_na), ~ !is.na(.)))
  # filter the data frame to keep only rows where all specified PRS columns are not NA
  
  out
}
#endregion

#region HELPERS
# ────────────────────────────────────────────────────────────────────────────────
# 5. Remove outliers based on Z-score ----

# this function removes outliers from numeric variables based on Z-score threshold
# useful for cleaning PRS data before analysis, controls for errors

remove_outliers_z <- function(df, numeric_vars, thresh = 3) {
    # takes a data frame, a vector of numeric variable names, and a Z-score threshold
  stopifnot(is.data.frame(df), is.numeric(thresh), thresh > 0)
  # check that df is a data frame, thresh is numeric, and thresh > 0
  for (v in numeric_vars) {
    if (!v %in% names(df)) stop("Variable not found: ", v) # if variable is not found in df stop execution
    z <- (df[[v]] - mean(df[[v]], na.rm = TRUE)) / sd(df[[v]], na.rm = TRUE) # compute Z-scores (for each variable, divide value by sd; NA values are ignored)
    df <- df[abs(z) <= thresh, , drop = FALSE] # remove outliers in absolute Z-score
  }
  # for loop: remove outliers from each numeric variable (main script uses PRSs)
  df
}

.ensure_factor01 <- function(x, control_label = "control", case_label = "case") {
  # accept 0/1 numeric or any two-level factor or character; order is control then case
  if (is.numeric(x)) {
    factor(x, levels = c(0,1), labels = c(control_label, case_label)) # coerce numeric 0/1 to labeled factor
  } else {
    x <- as.factor(x) # coerce to factor for consistent handling
    if (nlevels(x) != 2) stop("status_col must have exactly 2 levels.") # validate binary outcome
    levels(x) <- c(control_label, case_label) # rename levels to control and case keeping order
    x
  }
}

# infer which level of a 2-level factor corresponds to the positive "case" class.
# Heuristics: look for common labels; otherwise assume second level is the case.
.infer_case_level <- function(levels_vec) {
  # decide which factor level corresponds to the positive case class
  lv <- tolower(levels_vec) # normalize for robust pattern matching
  hit <- which(grepl("case|yes|pos|positive|1|migraine", lv)) # search common positive labels
  if (length(hit) >= 1) return(hit[1]) # use the first matching level if present
  if (length(lv) >= 2) return(2L) # otherwise assume the second level is case
  return(1L) # fallback to first level if only one exists
}

# Convert outcome to 0/1 numeric with 1 = case (based on .infer_case_level for factors)
.as_binary01 <- function(y) {
  # convert outcome to 0 or 1 numeric with 1 representing the case class
  if (is.numeric(y)) {
    return(as.numeric(y)) # already numeric, coerce to ensure integer-like
  }
  if (is.factor(y)) {
    lv <- levels(y) # get factor levels
    case_idx <- .infer_case_level(lv) # decide which level is case
    return(as.integer(y == lv[case_idx])) # encode case as 1 and control as 0
  }
  yf <- as.factor(y) # fallback: coerce character to factor
  if (nlevels(yf) != 2) stop("Outcome must be binary to convert to 0/1") # enforce binary outcome
  lv <- levels(yf) # get factor levels
  case_idx <- .infer_case_level(lv) # decide which level is case
  as.integer(yf == lv[case_idx]) # encode case as 1 and control as 0
}

# Predict probabilities for the case class from a fitted glm(binomial) model.
.predict_case_prob <- function(fit, newdata = NULL) {
  # get predicted probability for the case class from a fitted glm binomial model
  mf_fit <- stats::model.frame(fit) # extract model frame used in fitting
  y_fit  <- stats::model.response(mf_fit) # extract original response to infer coding
  nd <- if (is.null(newdata)) mf_fit else newdata # choose data to predict on
  preds <- stats::predict(fit, type = "response", newdata = nd) # obtain predicted probabilities
  if (is.numeric(y_fit)) {
    return(preds) # numeric 0 or 1 returns P(y==1) directly
  }
  if (!is.factor(y_fit) || nlevels(y_fit) != 2) # ensure binary factor when not numeric
    stop("Outcome must be a 2-level factor or numeric 0/1")
  lv <- levels(y_fit)
  case_idx <- .infer_case_level(lv) # locate which level is case
  first_is_case <- identical(lv[1], lv[case_idx]) # check if first level equals case
  if (first_is_case) preds else (1 - preds) # invert if model predicted P(first level)
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 6. Plot sex distribution ----

# this function creates bar plot of sex distribution by migraine status
# useful for demographic visualization in EDA

#region PLOT SEX DISTRIBUTION
plot_theme_readable <- function(sz = 16) {
  # define a minimal theme with larger, readable text for all plot elements
  theme_minimal(base_size = sz) + # base theme with configurable font size
    theme(
      plot.title = element_blank(), # remove title by default
      plot.subtitle = element_blank(),  # remove subtitle by default
      axis.title.x = element_text(size = sz), # set x/y axis title and x/y tick label size
      axis.title.y = element_text(size = sz),
      axis.text.x  = element_text(size = sz),
      axis.text.y  = element_text(size = sz), 
      legend.title = element_text(size = sz), 
      legend.text  = element_text(size = sz) # set legend text size
    )
}

plot_sex_distribution <- function(df, status_col, sex_col) {
  # show sex distribution within outcome groups as stacked proportions
  ggplot(df %>% mutate(.status = .ensure_factor01(.data[[status_col]]), # coerce outcome to control and case labels
                       .sex = as.factor(.data[[sex_col]])), # ensure sex is treated as factor
         aes(x = .status, fill = .sex)) + # map outcome on x and sex to fill
    geom_bar(position = "fill", width = 0.6) + # stack to proportions within each group
    scale_y_continuous(labels = scales::percent_format()) + # show y axis as percentages
    labs(x = "Status", y = "Proportion", fill = sex_col) + # set axis and legend labels
    plot_theme_readable() # apply readable theme
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 7. Plot chip distribution ----

# this function creates bar plot of genotyping array distribution by migraine status
# helps identify potential batch effects from different arrays

#region PLOT CHIP DISTRIBUTION
plot_chip_distribution <- function(df, status_col, chip_col) {
  # show genotyping array distribution within outcome groups as stacked proportions
  ggplot(df %>% mutate(.status = .ensure_factor01(.data[[status_col]]), # coerce outcome to control and case labels
                       .chip = as.factor(.data[[chip_col]])), # ensure chip is treated as factor
         aes(x = .status, fill = .chip)) + # map outcome on x and chip to fill
    geom_bar(position = "fill", width = 0.6) + # stack to proportions within each group
    scale_y_continuous(labels = scales::percent_format()) + # show y axis as percentages
    labs(x = "Status", y = "Proportion", fill = chip_col) + # set axis and legend labels
    plot_theme_readable() # apply readable theme
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 8. Choose two-sample test ----

# this function selects appropriate statistical test based on data characteristics
# helper function for statistical testing in EDA

#region CHOOSE 2SAMPLE TEST
choose_2sample_test <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)] # remove non-finite values from both samples
  n_x <- length(x); n_y <- length(y) # store sample sizes for checks and reporting
  if (n_x < 2 || n_y < 2) { # require at least two observations per group
    return(list(test = "insufficient data", p.value = NA_real_, fit = NULL,
                n_x = n_x, n_y = n_y)) # return placeholder when there is not enough data
  }

  # detect severe outliers via 3*IQR fences
  has_severe_outliers <- function(v) {
    v <- v[is.finite(v)] # drop non-finite values before computing fences
    if (length(v) < 4) return(FALSE) # need enough points for IQR to be meaningful
    qs  <- stats::quantile(v, c(0.25, 0.75), na.rm = TRUE) # compute first and third quartiles
    iqr <- diff(qs) # compute interquartile range
    if (!is.finite(iqr) || iqr == 0) return(FALSE) # skip if dispersion is degenerate
    lower <- qs[1] - 3 * iqr # lower 3*IQR fence
    upper <- qs[2] + 3 * iqr # upper 3*IQR fence
    any(v < lower | v > upper, na.rm = TRUE) # flag presence of severe outliers
  }

  if (has_severe_outliers(x) || has_severe_outliers(y)) { # prefer robust test when severe outliers exist
    fit <- tryCatch(stats::wilcox.test(x, y, exact = FALSE), error = function(e) NULL) # rank-sum test without exact calc
    if (is.null(fit)) { # handle rare failures in wilcox.test
      return(list(test = "Wilcoxon failed", p.value = NA_real_, fit = NULL,
                  n_x = n_x, n_y = n_y)) # return failure marker with sample sizes
    }
    return(list(test = "Wilcoxon rank-sum (outliers)", p.value = fit$p.value, fit = fit,
                n_x = n_x, n_y = n_y)) # return robust test result and metadata
  }

  # default: Welch's t-test (robust to unequal variances)
  fit <- stats::t.test(x, y, var.equal = FALSE) # Welch t-test defaults to unequal variances
  list(test = "Welch t-test", p.value = fit$p.value, fit = fit, n_x = n_x, n_y = n_y) # return parametric test result
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 9. Comprehensive preliminary EDA ----

# this function performs comprehensive exploratory data analysis with effect sizes
# generates summary statistics, plots, and saves results to file

#region PRELIMIN EDA
# AI transparency: this function was heavily reviewed, improved and debugged with the help of Copilot
prelimin_eda <- function(df, status_var, sex_var, age_var, chip_var, female_code,  # core inputs: df and key column names
                         outfile = "prelimin_eda_summary.txt", output_dir,       # output file and directory
                         dataset_label = NULL,                                    # optional label for headers
                         n_with_outliers = NULL) {                                # N before outlier removal if known
  
  # Check variable presence
  vars <- c(status_var, sex_var, age_var, chip_var)  # required columns for summaries and tests
  if (!all(vars %in% names(df))) {  # validate all required columns exist
    stop("Missing required columns: ", paste(vars[!vars %in% names(df)], collapse=", "))  # fail fast if any missing
  }
  
  # Check that sex variable has exactly 2 levels (ignoring NA)
  sex_levels <- unique(na.omit(df[[sex_var]]))  # derive observed sex codes (ignoring NAs)
  if (length(sex_levels) != 2) {  # ensure binary sex coding for comparisons
    stop(sprintf("Sex variable '%s' must have exactly 2 non-NA levels, but found: %s",
                 sex_var, paste(sex_levels, collapse = ", ")))  # informative error for mis-coded sex
  }
  
  # Standardize status to "Control"/"Case"
  status_label <- if (is.numeric(df[[status_var]])) {  # standardize status to Control/Case
    factor(df[[status_var]], levels = c(0,1), labels = c("Control","Case"))  # enforce order 0=Control,1=Case
  } else {
    x <- as.factor(df[[status_var]])  # coerce to factor for consistent handling
    if (nlevels(x) != 2) stop("status_var must have exactly 2 levels.")  # require binary outcome
    levels(x) <- c("Control","Case")  # set readable labels
    x  # return labeled factor
  }
  
  # Coerce female_code to appropriate type
  sex_vec <- df[[sex_var]]  # raw sex vector as present in df
  female_code_coerced <- if (is.numeric(sex_vec)) suppressWarnings(as.numeric(female_code)) else as.character(female_code)  # align type of code
  message(sprintf("Using '%s' as female code for variable '%s'.", female_code_coerced, sex_var))  # explicit log for reproducibility
  
  # Build a clean sex label using coerced female code
  sex_label <- ifelse(sex_vec == female_code_coerced, "Female", "Male")  # map to Male/Female labels
  sex_label <- factor(sex_label, levels = c("Male","Female"))  # fix ordering: Male first
  
  outfile <- file.path(output_dir, outfile)  # build absolute output path
  
  df <- dplyr::mutate(df,  # augment df with standardized helper columns used below
                      .status = status_label,  # "Control"/"Case" factor
                      .sex    = sex_label,     # "Male"/"Female" factor
                      .age    = .data[[age_var]],  # numeric age for summaries/tests
                      .chip   = as.factor(.data[[chip_var]])  # categorical genotyping array
  )
  
  df$sex_label <- factor(ifelse(df[[sex_var]] == female_code_coerced, "Female", "Male"),  # legacy column name used in downstream code
                         levels = c("Male","Female"))
  
  # Filtered subsets
  df_age            <- df[!is.na(df[[age_var]]), ]  # rows with age present
  df_age_sex        <- df[!is.na(df[[age_var]]) & !is.na(df[[sex_var]]), ]  # age + sex present
  df_age_status     <- df[!is.na(df[[age_var]]) & !is.na(df[[status_var]]), ]  # age + status present
  df_age_sex_status <- df[!is.na(df[[age_var]]) & !is.na(df[[sex_var]]) & !is.na(df[[status_var]]), ]  # age + sex + status present
  
  n_total    <- nrow(df)  # total N used in summaries
  pct_female <- mean(df$sex_label == "Female", na.rm = TRUE) * 100  # % female
  mean_age   <- mean(df_age[[age_var]])  # mean age (non-missing)
  sd_age     <- sd(df_age[[age_var]])  # sd age
  min_age    <- min(df_age[[age_var]])  # min age
  max_age    <- max(df_age[[age_var]])  # max age
  
  get_stat_summary <- function(data, condition) {  # helper: min/max/mean under a logical mask
    vals <- data[[age_var]][condition]  # extract values for condition
    if (length(vals) == 0 || all(is.na(vals))) return(c(NA, NA, NA))  # guard small/NA-only slices
    c(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE), mean(vals, na.rm = TRUE))  # return trio
  }
  
  # Get summary statistics for females
  female_stats <- get_stat_summary(df_age_sex, df_age_sex[[sex_var]] == female_code_coerced)  # stats among females
  min_age_fem <- female_stats[1]
  max_age_fem <- female_stats[2]
  mean_age_fem <- female_stats[3]
  
  # Get summary statistics for males
  male_stats <- get_stat_summary(df_age_sex, df_age_sex[[sex_var]] != female_code_coerced)  # stats among males
  min_age_male <- male_stats[1]
  max_age_male <- male_stats[2]
  mean_age_male <- male_stats[3]
  
  # Use the unified chooser: Welch's t-test by default; Wilcoxon if severe outliers
  safe_choose <- function(x, y) {  # robust wrapper for test selection with NA/size guards
    if (length(na.omit(x)) < 2 || length(na.omit(y)) < 2) {  # need >=2 per group
      return(list(fit = list(p.value = NA_real_), test = "Not enough data"))  # mark insufficient data
    }
    tryCatch(choose_2sample_test(x, y),  # prefer Welch; fallback to Wilcoxon under severe outliers
             error = function(e) list(fit = list(p.value = NA_real_), test = "Error"))  # swallow rare errors
  }

  res_sex_overall   <- safe_choose(df_age_sex[[age_var]][df_age_sex[[sex_var]] == female_code_coerced],  # F vs M overall
                                   df_age_sex[[age_var]][df_age_sex[[sex_var]] != female_code_coerced])
  res_cases         <- safe_choose(df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] == female_code_coerced & df_age_sex_status$.status == "Case"],
                                   df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] != female_code_coerced & df_age_sex_status$.status == "Case"])
  res_controls      <- safe_choose(df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] == female_code_coerced & df_age_sex_status$.status == "Control"],
                                   df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] != female_code_coerced & df_age_sex_status$.status == "Control"])
  res_age_cs        <- safe_choose(df_age_status[[age_var]][df_age_status$.status == "Case"],  # Case vs Control overall
                                   df_age_status[[age_var]][df_age_status$.status == "Control"])
  res_m             <- safe_choose(df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Male"   & df_age_sex_status$.status == "Case"],
                                   df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Male"   & df_age_sex_status$.status == "Control"])
  res_f             <- safe_choose(df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Female" & df_age_sex_status$.status == "Case"],
                                   df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Female" & df_age_sex_status$.status == "Control"]) 
  
  sex_tab <- table(df_age_sex_status$.status, df_age_sex_status$sex_label)  # 2x2 contingency: status by sex
  if (nrow(sex_tab) < 2 || ncol(sex_tab) < 2) {  # guard degenerate tables
    # if < 2 rows or columns, not enough data for test
    sex_p <- NA_real_  # no p-value
    sex_test <- "Insufficient dimensions"  # note lack of dimensions
    cramers_v_sex <- NA_real_  # effect size not defined
  } else if (any(sex_tab < 5)) {  # small counts → Fisher's exact
    sex_p <- fisher.test(sex_tab)$p.value  # exact test p-value
    sex_test <- "Fisher's exact"  # record test name
    # Cramér's V calculation
    chi2_sex <- chisq.test(sex_tab)$statistic  # use chi-square statistic for V
    n_sex <- sum(sex_tab)  # total N
    cramers_v_sex <- sqrt(chi2_sex / (n_sex * (min(dim(sex_tab)) - 1)))  # effect size
  } else {
    chi2_test_sex <- chisq.test(sex_tab)  # standard chi-squared test
    sex_p <- chi2_test_sex$p.value  # p-value from test
    sex_test <- "Chi-squared"  # record test name
    # Cramér's V calculation
    n_sex <- sum(sex_tab)  # total N
    cramers_v_sex <- sqrt(chi2_test_sex$statistic / (n_sex * (min(dim(sex_tab)) - 1)))  # effect size
  }
  
  chip_tab <- table(df_age_sex_status$.status, df_age_sex_status[[chip_var]])  # status by genotyping chip
  if (nrow(chip_tab) < 2 || ncol(chip_tab) < 2) {  # guard degenerate tables
    chip_p <- NA_real_  # no p-value
    chip_test <- "Insufficient dimensions"  # record limitation
    cramers_v_chip <- NA_real_  # effect size not defined
  } else if (any(chip_tab < 5)) {  # small expected counts → Fisher's exact
    chip_p <- fisher.test(chip_tab)$p.value  # exact p-value
    chip_test <- "Fisher's exact"  # record method
    # Cramér's V calculation for Fisher's test
    chi2_chip <- chisq.test(chip_tab)$statistic  # use chi-square stat for V
    n_chip <- sum(chip_tab)  # total N
    cramers_v_chip <- sqrt(chi2_chip / (n_chip * (min(dim(chip_tab)) - 1)))  # effect size
  } else {
    chi2_test_chip <- chisq.test(chip_tab)  # standard chi-squared
    chip_p <- chi2_test_chip$p.value  # p-value
    chip_test <- "Chi-squared"  # record method
    # Cramér's V calculation
    n_chip <- sum(chip_tab)  # total N
    cramers_v_chip <- sqrt(chi2_test_chip$statistic / (n_chip * (min(dim(chip_tab)) - 1)))  # effect size
  }
  
  # Calculate descriptive statistics for better interpretation
  age_overall_fem <- mean(df_age_sex[[age_var]][df_age_sex[[sex_var]] == female_code_coerced], na.rm = TRUE)  # mean age females
  age_overall_mal <- mean(df_age_sex[[age_var]][df_age_sex[[sex_var]] != female_code_coerced], na.rm = TRUE)  # mean age males
  age_diff_sex <- age_overall_fem - age_overall_mal  # difference females - males
  
  age_cases_overall <- mean(df_age_status[[age_var]][df_age_status$.status == "Case"], na.rm = TRUE)  # mean age cases
  age_controls_overall <- mean(df_age_status[[age_var]][df_age_status$.status == "Control"], na.rm = TRUE)  # mean age controls
  age_diff_status <- age_cases_overall - age_controls_overall  # difference cases - controls
  
  # Additional comparisons
  age_fem_cases <- mean(df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] == female_code_coerced & df_age_sex_status$.status == "Case"], na.rm = TRUE)  # mean age F cases
  age_mal_cases <- mean(df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] != female_code_coerced & df_age_sex_status$.status == "Case"], na.rm = TRUE)  # mean age M cases
  age_diff_sex_cases <- age_fem_cases - age_mal_cases  # difference F-M among cases
  
  age_fem_controls <- mean(df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] == female_code_coerced & df_age_sex_status$.status == "Control"], na.rm = TRUE)  # mean age F controls
  age_mal_controls <- mean(df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] != female_code_coerced & df_age_sex_status$.status == "Control"], na.rm = TRUE)  # mean age M controls
  age_diff_sex_controls <- age_fem_controls - age_mal_controls  # difference F-M among controls
  
  age_mal_cases_vs_controls <- mean(df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Male" & df_age_sex_status$.status == "Case"], na.rm = TRUE)  # mean age M cases
  age_mal_controls_vs_cases <- mean(df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Male" & df_age_sex_status$.status == "Control"], na.rm = TRUE)  # mean age M controls
  age_diff_mal_status <- age_mal_cases_vs_controls - age_mal_controls_vs_cases  # difference cases - controls in males
  
  age_fem_cases_vs_controls <- mean(df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Female" & df_age_sex_status$.status == "Case"], na.rm = TRUE)  # mean age F cases
  age_fem_controls_vs_cases <- mean(df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Female" & df_age_sex_status$.status == "Control"], na.rm = TRUE)  # mean age F controls
  age_diff_fem_status <- age_fem_cases_vs_controls - age_fem_controls_vs_cases  # difference cases - controls in females
  
  # Helper function for Cohen's d calculation
  calc_cohens_d <- function(x1, x2) {  # standardized mean difference (Hedges' g approx for large N)
    n1 <- length(x1[!is.na(x1)])  # sample size group 1
    n2 <- length(x2[!is.na(x2)])  # sample size group 2
    if (n1 < 2 || n2 < 2) return(NA)  # require at least two per group
    
    m1 <- mean(x1, na.rm = TRUE)  # group 1 mean
    m2 <- mean(x2, na.rm = TRUE)  # group 2 mean
    s1 <- var(x1, na.rm = TRUE)   # group 1 variance
    s2 <- var(x2, na.rm = TRUE)   # group 2 variance
    
    pooled_sd <- sqrt(((n1-1)*s1 + (n2-1)*s2) / (n1 + n2 - 2))  # pooled sd
    (m1 - m2) / pooled_sd  # Cohen's d
  }
  
  # Calculate Cohen's d for all comparisons
  cohens_d_sex <- calc_cohens_d(df_age_sex[[age_var]][df_age_sex[[sex_var]] == female_code_coerced],  # d for F vs M overall
                                df_age_sex[[age_var]][df_age_sex[[sex_var]] != female_code_coerced])
  
  cohens_d_status <- calc_cohens_d(df_age_status[[age_var]][df_age_status$.status == "Case"],  # d for Case vs Control overall
                                   df_age_status[[age_var]][df_age_status$.status == "Control"])
  
  cohens_d_sex_cases <- calc_cohens_d(df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] == female_code_coerced & df_age_sex_status$.status == "Case"],  # d for F vs M among cases
                                      df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] != female_code_coerced & df_age_sex_status$.status == "Case"])
  
  cohens_d_sex_controls <- calc_cohens_d(df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] == female_code_coerced & df_age_sex_status$.status == "Control"],  # d for F vs M among controls
                                         df_age_sex_status[[age_var]][df_age_sex_status[[sex_var]] != female_code_coerced & df_age_sex_status$.status == "Control"])
  
  cohens_d_mal_status <- calc_cohens_d(df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Male" & df_age_sex_status$.status == "Case"],  # d for Case vs Control in males
                                       df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Male" & df_age_sex_status$.status == "Control"])
  
  cohens_d_fem_status <- calc_cohens_d(df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Female" & df_age_sex_status$.status == "Case"],  # d for Case vs Control in females
                                       df_age_sex_status[[age_var]][df_age_sex_status$sex_label == "Female" & df_age_sex_status$.status == "Control"])
  
  # Header with dataset information and Ns
  if (!is.null(n_with_outliers)) {  # include outlier removal counts when provided
    removed_n <- n_with_outliers - n_total # compute number of removed outliers
    header <- c( # this creates the header vector with outlier info
      sprintf("Dataset used for EDA: %s", if (is.null(dataset_label)) "unspecified" else dataset_label), # if dataset_label is NULL write "unspecified"
      sprintf("N (with outliers): %d", n_with_outliers), 
      sprintf("N (without outliers): %d", n_total),
      sprintf("Outliers removed: %d", removed_n),
      ""
    )
  } else {
    header <- c(  # simpler header when no outlier count is available
      sprintf("Dataset used for EDA: %s", if (is.null(dataset_label)) "unspecified" else dataset_label), # if dataset_label is NULL write "unspecified"
      sprintf("N (dataset used): %d", n_total),
      ""
    )
  }

    # header is a character vector containing dataset information and sample sizes, lines is a character vector containing summary statistics

  lines <- c(header, # builds character vector with formatted summary statistics via sprintf with the header at the top
    sprintf("N total: %d", n_total),
    sprintf("%% female: %.1f", pct_female),
    sprintf("Mean age: %.2f", mean_age),
    sprintf("SD age: %.2f", sd_age),
    sprintf("Age range: %.1f - %.1f", min_age, max_age),
    sprintf("Age range females: %.1f - %.1f", min_age_fem, max_age_fem),
    sprintf("Age range males: %.1f - %.1f", min_age_male, max_age_male),
    "", # this adds a blank line for separation
    "=== AGE COMPARISONS ===",
    sprintf("Age F vs M (overall): F=%.2f vs M=%.2f (diff=%.2f yrs), p=%.2e, Cohen's d=%.3f, test=%s", 
            age_overall_fem, age_overall_mal, age_diff_sex, res_sex_overall$fit$p.value, cohens_d_sex, res_sex_overall$test),
    sprintf("Age Case vs Control (overall): Cases=%.2f vs Controls=%.2f (diff=%.2f yrs), p=%.2e, Cohen's d=%.3f, test=%s", 
            age_cases_overall, age_controls_overall, age_diff_status, res_age_cs$fit$p.value, cohens_d_status, res_age_cs$test),
    sprintf("Age F vs M (cases): F=%.2f vs M=%.2f (diff=%.2f yrs), p=%.2e, Cohen's d=%.3f, test=%s", 
            age_fem_cases, age_mal_cases, age_diff_sex_cases, res_cases$fit$p.value, cohens_d_sex_cases, res_cases$test),
    sprintf("Age F vs M (controls): F=%.2f vs M=%.2f (diff=%.2f yrs), p=%.2e, Cohen's d=%.3f, test=%s", 
            age_fem_controls, age_mal_controls, age_diff_sex_controls, res_controls$fit$p.value, cohens_d_sex_controls, res_controls$test),
    sprintf("Age Case vs Control in males: Cases=%.2f vs Controls=%.2f (diff=%.2f yrs), p=%.2e, Cohen's d=%.3f, test=%s", 
            age_mal_cases_vs_controls, age_mal_controls_vs_cases, age_diff_mal_status, res_m$fit$p.value, cohens_d_mal_status, res_m$test),
    sprintf("Age Case vs Control in females: Cases=%.2f vs Controls=%.2f (diff=%.2f yrs), p=%.2e, Cohen's d=%.3f, test=%s", 
            age_fem_cases_vs_controls, age_fem_controls_vs_cases, age_diff_fem_status, res_f$fit$p.value, cohens_d_fem_status, res_f$test),
    "",
    "=== CATEGORICAL COMPARISONS ===",
    sprintf("Sex distribution by status: p = %.2e, Cramér's V = %.3f, test = %s", sex_p, cramers_v_sex, sex_test),
    sprintf("   -> Interpretation: Tests whether migraine cases and controls differ in sex distribution"),
    sprintf("Chip distribution by status: p = %.3g, Cramér's V = %.3f, test = %s", chip_p, cramers_v_chip, chip_test),
    sprintf("   -> Interpretation: Tests whether migraine cases and controls differ across genotyping chips"),
    "",
    "=== INTERPRETATION NOTES ===",
    "Cohen's d interpretation: <0.2=negligible, 0.2=small, 0.5=medium, 0.8=large effect",
    "Cramér's V interpretation: <0.1=negligible, 0.1=small, 0.3=medium, 0.5=large association",
    "For categorical tests: p<0.05 = statistically significant, BUT focus on effect sizes",
    "With N>140k, tiny differences become statistically significant",
    "Focus on effect sizes and practical significance, not just p-values",
    "Consider your research context to determine what differences are meaningful",
    "More info: https://statisticsbyjim.com/hypothesis-testing/practical-statistical-significance/"
  )
  
  # write or append the summary
  if (file.exists(outfile)) {  # append with a separator if file pre-exists
    cat("\n", paste(rep("-", 80), collapse = ""), "\n", file = outfile, append = TRUE)  # this line adds a visual delimiter which is a series of dashes
    cat(paste0(paste(lines, collapse = "\n"), "\n"), file = outfile, append = TRUE)  # append block
  } else {
    writeLines(lines, con = outfile)  # write new file
  }
  
  # in main script, results are written for without_outliers first and the with_outliers is appended
  
  invisible(NULL)  # no return needed
}

# ────────────────────────────────────────────────────────────────────────────────
# 10. Plot PRS decile prevalence ----

# this function creates decile-prevalence plot for PRS variable
# shows risk gradient across PRS distribution

#region PLOT PRS DECILE PREVALENCE
plot_prs_decile <- function(df, prs_var, status_var) {
  df %>%
    mutate(decile = ntile(.data[[prs_var]], 10)) %>% # calculates deciles for the column in the current data frame whose name is in prs_var 
                                                     # by dividing the data into 10 equal parts
    group_by(decile) %>% # group by decile
    summarise( #summarise() creates a new dataframe
      prev = mean(.data[[status_var]]), # calculates mean prevalence of outcome variable for each decile
      se   = sqrt(prev * (1 - prev) / n()) # calculates standard error for each decile
    ) %>%
  ggplot(aes(decile, prev)) + 
    geom_point() +
    geom_errorbar(aes(ymin = prev - 1.96 * se, # calculates 95% CI for each decile; the use of 1.96 is explained in the report
                      ymax = prev + 1.96 * se), width = .2) +
  labs(x = "PRS decile", y = "Migraine prevalence") +
  plot_theme_readable()
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 11. Plot univariate odds ratios ----

# this function creates forest plot of univariate odds ratios for covariates
# useful for covariate screening and effect size visualization

#region PLOT UNIVARIABLE OR
plot_uni_or <- function(df, covars, status_var, outfile = NULL) {

  .esc <- function(s) gsub("([][{}()+*^$|\\.?\\\\])", "\\\\\\\\\\\\1", s) # escape special characters for regex as a precaution

  ors <- purrr::map_dfr(covars, function(v) {
    if (!v %in% names(df)) return(NULL) # if variable is not in dataframe return NULL
    x <- df[[v]] # extract variable from dataframe
    fit <- tryCatch( # fit logistic regression model with tryCatch for error handling
      glm(as.formula(sprintf("%s ~ `%s`", status_var, v)), data = df, family = binomial()), # sprintf() is used to create the formula with
                                                                                # outcome variable and variable names; %s is a placeholder for variable names
      error = function(e) NULL # if model fitting fails, return NULL
    )
    if (is.null(fit)) return(NULL) # if model is NULL, return NULL

    tt <- broom::tidy(fit) %>% filter(term != "(Intercept)") # tidy the model output and keep term that is not intercept
    if (nrow(tt) == 0) return(NULL) # if tidy output is empty, return NULL

    is_cat <- is.factor(x) || is.character(x) # check if variable is categorical
    ref <- if (is_cat) levels(as.factor(x))[1] else NA_character_ # get reference level for categorical variables
    lvl <- if (is_cat) { # get levels for categorical variables
      gsub(paste0("^`?", .esc(v), "`?"), "", tt$term) # remove variable name from term
    } else {
      rep("", nrow(tt)) # else, return empty string for each level, this is for non-categorical variables
    }
    label <- if (is_cat) paste0(v, ": ", lvl, " vs ", ref) else v # create label for plot using variable name and levels for categorical variables

    tt %>% 
      mutate(
        variable = v,
        label = label,
        OR   = exp(estimate),
        LCL  = exp(estimate - 1.96 * std.error),
        UCL  = exp(estimate + 1.96 * std.error)
      ) %>%
      select(variable, label, estimate, std.error, statistic, p.value, OR, LCL, UCL)
  })

  if (is.null(ors) || nrow(ors) == 0) { # if no valid ORs were computed raise a warning and return NULL
    warning("No valid univariable ORs were computed.")
    return(invisible(NULL))
  }

  p <- ggplot(ors, aes(x = label, y = OR, ymin = LCL, ymax = UCL)) +
    geom_pointrange() + geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() + plot_theme_readable() +
  labs(x = "", y = "Unadjusted OR (95% CI)")
  if (!is.null(outfile)) ggsave(outfile, p, width = 7, height = 4, dpi = 300)
  p
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 12. Save PRS decile plots ----

# this function saves PRS decile-prevalence plot to file
# useful for visualizing PRS-outcome relationships across risk spectrum

#region EDA SAVE PRS DECILE
eda_save_prs_decile <- function(df, prs_var, status_var, outdir,
                                width = 6, height = 4, dpi = 300) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  p <- plot_prs_decile(df, prs_var, status_var)
  ggsave(file.path(outdir, paste0("decile_prev_", prs_var, ".png")),
         p, width = width, height = height, dpi = dpi)
  invisible(p)
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 13. Save univariate OR plots ----

# this function saves univariate OR forest plot to file
# creates standardized forest plots for covariate associations

#region EDA SAVE UNI OR
eda_save_uni_or <- function(df, covariates, status_var, outdir,
                            width = 7, height = 4, dpi = 300) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE) # if output directory doesn't exist, create it
  p <- plot_uni_or(df, covariates, status_var, outfile = NULL) # create univariate OR plot
  ggsave(file.path(outdir, "univ_or.png"), p, width = width, height = height, dpi = dpi) # save plot to file
  invisible(p)
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 14. Run all EDA functions ----

# this function runs complete EDA pipeline and saves all plots
# comprehensive wrapper for all exploratory analysis functions

#region RUN ALL EDA

run_all_eda <- function(ukbb_clean, prs_cols, status_col, covariates,  
                        output_dir,    
                        facet_vars, # optional: names for sex and chip columns
                        subdir = "EDA", # subfolder for EDA outputs
                        viz_z_thresh = 5) { # visualization-only Z cutoff
  # threshold of 5 used to be conservative
  stopifnot(is.data.frame(ukbb_clean))  # ensure input is a data.frame
  if (length(prs_cols) < 2) stop("Expected at least two PRS columns in prs_cols.")  # need two PRS for plots
  if (!status_col %in% names(ukbb_clean)) stop("status_col not in data.")  # status must exist
  if (!all(prs_cols %in% names(ukbb_clean))) stop("Some prs_cols not in data.")  # all PRS columns must exist
  
  outdir <- file.path(output_dir, subdir)  # EDA output directory path
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)  # create if missing
  
  p1 <- prs_cols[1]  # first PRS for plotting
  p2 <- prs_cols[2]  # second PRS for plotting
  out_files <- list()  # collect output file paths for downstream use
  
  df_plot <- ukbb_clean  # copy to avoid modifying caller's data

  # Drop duplicate outlier versions to avoid confusion (main dataframe only)
  
  # decile plots: prevalence by PRS decile for both PRS
  eda_save_prs_decile(ukbb_clean, p1, status_col, outdir)  # save decile plot for first PRS
  eda_save_prs_decile(ukbb_clean, p2, status_col, outdir)  # save decile plot for second PRS
  out_files$decile <- file.path(outdir, c(paste0("decile_prev_", p1, ".png"),  # record paths
                                          paste0("decile_prev_", p2, ".png")))
  
  # univariate OR plot: forest plot for covariates
  eda_save_uni_or(ukbb_clean, covariates, status_col, outdir)  # save univariable OR plot
  out_files$uni_or <- file.path(outdir, "univ_or.png")  # record path
  
  # Sex and chip categorical distributions by status (optional if vars provided)
  sex_fv  <- if (length(facet_vars) >= 1) facet_vars[[1]] else NULL  # pick sex column name
  chip_fv <- if (length(facet_vars) >= 2) facet_vars[[2]] else NULL  # pick chip column name

  if (!is.null(sex_fv) && sex_fv %in% names(ukbb_clean)) {  # guard column existence
    p_sex <- plot_sex_distribution(ukbb_clean, status_col, sex_fv)  # build plot
    ggsave(file.path(outdir, "sex_by_status.png"), p_sex, width = 10, height = 6, dpi = 300)  # save
    out_files$sex_by_status <- file.path(outdir, "sex_by_status.png")  # record path
  }
  if (!is.null(chip_fv) && chip_fv %in% names(ukbb_clean)) {  # guard column existence
    p_chip <- plot_chip_distribution(ukbb_clean, status_col, chip_fv)  # build plot
    ggsave(file.path(outdir, "chip_by_status.png"), p_chip, width = 10, height = 6, dpi = 300)  # save
    out_files$chip_by_status <- file.path(outdir, "chip_by_status.png")  # record path
  }
  
  message("EDA plots saved to: ", outdir)  # notify user where outputs are saved
  invisible(out_files)  # return list of output paths invisibly
}

# ────────────────────────────────────────────────────────────────────────────────
# 15. Fit PRS association model ----

# this function fits a logistic regression model to assess the association
# between two PRS variables and a binary outcome, optionally including covariates;
# it allows for an interaction term between the two PRS variables
# and returns a tidy summary of the results, calculating odds ratios for easier
# interpretation and confidence intervals; controls for errors

#region FIT PRS MODEL
# AI transparency: this function was heavily reviewed, improved and debugged with the help of Copilot
fit_prs_model <- function(df, outcome_var, prs_vars,
                          covariates=character(), interaction=FALSE, 
                          sex_var = NULL, 
                          sex_interaction = FALSE,
                          alpha = 0.05) {
  
  # basic sanity checks
  all_vars <- c(outcome_var, prs_vars, covariates)
  missing <- setdiff(all_vars, names(df))
  if (length(missing)) stop("Columns not in df: ", paste(missing, collapse = ", "))
  # check if all specified variables are present in the data frame
  
  # outcome must be 0/1 or two-level factor
  y <- df[[outcome_var]]
  if (!(is.numeric(y) && all(unique(y) %in% c(0,1))) &&
      !(is.factor(y) && nlevels(y) == 2)) {
    stop("`", outcome_var, "` must be numeric 0/1 or a two‐level factor")
  }
  # if the outcome variable is not numeric 0/1 or a two-level factor, stop execution
  
  message(">>> Running association analysis …")
  
  terms <- prs_vars # this is a vector of PRS variable names
  
  # if PRSxPRS term
  if (interaction) { # if interaction is TRUE, add interaction term
    terms <- c(terms, paste0(prs_vars[1],":",prs_vars[2]))
  }
  
  # if PRSxsex term
  if (sex_interaction) {
    if (is.null(sex_var) || ! sex_var %in% names(df)) {
      stop("You asked for sex_interaction but ",
           sex_var, " is missing in df")
    }
    # add main effect of sex (if not already in covariates)
    if (! sex_var %in% covariates) {
      covariates <- c(covariates, sex_var)
    }
    # add each PRS×sex
    sex_int_terms <- paste0(prs_vars, ":", sex_var)
    terms <- c(terms, sex_int_terms)
  }
  
  terms <- c(terms, covariates)
  
  formstr <- sprintf("%s ~ %s", # create formula string using terms and covariates
                     outcome_var,
                     paste(terms, collapse = " + "))
  message("Fitting: ", formstr)
  
  form <- as.formula(formstr) # convert string to formula object; converting the string into a formula object allows us to use it in the glm() function
  # while leaving room for flexibility in specifying the model formula if needed 
  
  model <- glm(form, data = df, family = "binomial", na.action = na.exclude) # fit logistic regression model with NA handling
  
  tbl   <- broom::tidy(model) # tidy the model output
  
  #  inter is the character vector of coefficient names that you want to pull out of the model results
  inter <- prs_vars
  
  if (interaction) { # this adds PRS×PRS if requested
    inter <- c(inter, paste0(prs_vars[1], ":", prs_vars[2]))
  } 
  
  # add PRS×Sex if requested
  if (sex_interaction) {
    # this creates c("prs_asd_std:sex_plink", "prs_adhd_std:sex_plink")
    sex_ints <- paste0(prs_vars, ":", sex_var)
    inter     <- c(inter, sex_ints)
  }
  
  res <- tbl %>%
    dplyr::filter(term %in% inter) %>%
    dplyr::mutate(
      OR        = exp(estimate),
      # Logistic regression scale vs. OR: in a logistic regression, the estimate is 
      # on the log-odds scale: a 1-unit increase in PRS changes the log-odds of 
      # migraine by β. exp() converts that to an odds ratio, which is more interpretable
      
      LCL       = exp(estimate - 1.96 * std.error),
      # the 95% CI is computed on the log-odds scale as β±1.96×SE and then exponentiate 
      # those limits to get a CI for the OR; we use 1.96 because Z=±1.96 are the 2.5th 
      # and 97.5th percentiles of N(0,1); A 95% confidence interval for the standard 
      # normal distribution, then, is the interval (-1.96, 1.96), since 95% of the area 
      # under the curve falls within this interval 
      # (https://math.stackexchange.com/questions/1480904/given-a-95-confidence-interval-why-are-we-using-1-96-and-not-1-64)
      
      UCL       = exp(estimate + 1.96 * std.error),
      p.adj     = p.adjust(p.value, method = "bonferroni"), # Bonferroni correction for multiple testing
      threshold = alpha / length(inter), # Bonferroni correction threshold
      significant = p.adj < threshold # determine significance based on adjusted p-value
    ) %>%
    dplyr::select(term, estimate, std.error, p.value, p.adj, threshold, significant, OR, LCL, UCL)
  
  sig_terms <- res$term[res$significant]
  # extract significant terms based on adjusted p-value
  
  message("🎉 Association analysis complete.")
  
  list(
    model      = model,
    results    = res,
    pvalues    = setNames(res$p.value, res$term),
    p.adj      = setNames(res$p.adj, res$term),
    threshold  = unique(res$threshold),
    sig_terms  = sig_terms
  ) # return a list containing the model, results, p-values, adjusted p-values, 
  # threshold, and significant terms
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 15b. Run PRS model stratified by sex (no interactions; sex not in covariates) ----

# This helper runs the same logistic model separately within each sex level.
# It removes the sex variable from covariates (if present) and forces
# interaction = FALSE and sex_interaction = FALSE. It also skips strata with
# no variation in the outcome.

#region RUN SEX STRATIFIED (NO INTERACTIONS)
run_sex_stratified_no_interactions <- function(df, 
                                               outcome_var,
                                               prs_vars, 
                                               covariates = character(),# other covariates (sex removed per stratum)
                                               sex_var) { # column indicating sex strata
  stopifnot(is.data.frame(df))  # type check for input
  needed <- c(outcome_var, prs_vars, sex_var)  # required columns per strata fit
  miss <- setdiff(needed, names(df))  # identify any missing columns
  if (length(miss)) stop("Columns not in df: ", paste(miss, collapse = ", "))  # fail fast if missing
  if (length(prs_vars) < 1) stop("Need at least one PRS variable in prs_vars")  # require at least one PRS

  sx <- df[[sex_var]]  # vector of sex values
  lvls <- unique(stats::na.omit(as.vector(sx)))  # observed non-NA strata levels
  if (length(lvls) < 2) warning(sprintf("Sex variable '%s' has <2 non-NA levels.", sex_var))  # warn on degenerate strata

  covars_no_sex <- setdiff(covariates, sex_var)  # ensure sex not included as covariate within-stratum

  res <- list()  # collect fits by level name
  for (lv in lvls) {  # loop over strata
    dsub <- df[!is.na(df[[sex_var]]) & df[[sex_var]] == lv, , drop = FALSE]  # subset to current sex level
    y <- dsub[[outcome_var]]  # outcome vector in stratum
    y_no_na <- stats::na.omit(y)  # remove missing
    if (length(unique(y_no_na)) < 2) {  # skip if no variation in outcome
      warning(sprintf("Skipping stratum '%s' (no outcome variation)", as.character(lv)))
      res[[as.character(lv)]] <- NULL  # record as NULL
      next  # move to next stratum
    }
    fit <- tryCatch(
      fit_prs_model(
        df = dsub,  # stratum-specific data
        outcome_var = outcome_var, # same outcome
        prs_vars = prs_vars,# same PRS predictors
        covariates = covars_no_sex, # covariates excluding sex
        interaction = FALSE, # no PRS×PRS interaction
        sex_var = NULL, # no sex variable in model
        sex_interaction = FALSE # no PRS×sex interactions
      ),
      error = function(e) {  # handle failures per stratum without aborting loop
        warning(sprintf("Stratum '%s' failed: %s", as.character(lv), conditionMessage(e)))
        NULL
      }
    )
    res[[as.character(lv)]] <- fit  # store result (list from fit_prs_model or NULL)
  }
  invisible(res)  # return named list of fits invisibly
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 16. Assumptions

# 16a. Component + residual plots ----

# this function draws component+residual plot for one predictor in a glm fit
# useful for checking linearity assumption in logistic regression

#region PLOT COMPONENT + RESIDUAL
plot_cr_gg <- function(fit, var) {
  beta <- coef(fit)[var] # extracts coefficient for the predictor
  x <- fit$model[[var]] # extracts predictor values
  res <- residuals(fit, type = "deviance") # extracts deviance residuals

  # align lengths in case residuals were padded by na.exclude (this means that <NA> values were removed)
  if (length(res) != length(x)) { # if lengths differ, keep non-NA
    res <- res[!is.na(res)]
  }
  n <- min(length(x), length(res)) # determine minimum length to avoid index issues
  df <- tibble( # constructs a tibble
    x     = x[seq_len(n)], # variable x is assigned the first n elements of an object x
    resid = res[seq_len(n)] # column resid is similarly assigned the first n elements of the object res
  ) %>% mutate(partial = resid + beta * x) # partial column is computed as the sum of resid and beta * x
  
  ggplot(df, aes(x = x, y = partial)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      x     = var,
      y     = "Partial Residual"
    ) +
    plot_theme_readable()
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 16b. VIF check for multicollinearity ----

#region VIF CHECK
# AI transparency: this function was heavily reviewed, improved and debugged with the help of Copilot
vif_check <- function(fit, thresh = 5, type = c("terms", "predictor")) {
  type <- match.arg(type) # match the argument to one of the allowed types

  if (inherits(fit, "glm")) { # if the model is a generalized linear model
    X <- model.matrix(fit) # create model matrix
    # drop intercept term if present
    if ("(Intercept)" %in% colnames(X)) {
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE] # colnames() retrieves column names of X (matrix)
    }
    
    # formula for VIF is 1 / (1 - r²), where r² is the coefficient of determination (R-squared) from the regression

    vifs <- sapply(colnames(X), function(var) { # sapply() applies a function to each column name and creates a vector of results
      other <- setdiff(colnames(X), var) # selects all columns of X except the current variable, representing independent variables in the model
      r2 <- summary( # computes regression summary, and $r.squared extracts the R^2 value of the model
        lm(X[, var] ~ X[, other], data = as.data.frame(X)) # fits a linear model where var is the dependent variable and the other are the independent variables
      )$r.squared
      1 / (1 - r2) # uses formula ( \text{VIF} = \frac{1}{1 - R^2} ) to compute the VIF for each predictor
    })
    # results in a named vector, where each entry corresponds to the VIF of a variable in X

  } else {
    # fallback for lm or other supported models
    vifs <- car::vif(fit, type = type)
  }
  
  list(
    vif_values = vifs,
    pass       = all(vifs < thresh)
  )
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 16c. Cook’s distance to test influential data points ----

#region COOKS DISTANCE CHECK
cooks_check <- function(fit, cutoff = NULL) {
  d <- cooks.distance(fit)
  # n <- length(d)
  # k <- length(coef(fit)) # number of predictors
  
  # if user does not supply a custom cut‑off, use 1
  if (is.null(cutoff)) cutoff <- 1
  # https://stats.stackexchange.com/questions/87962/cooks-distance-cut-off-value
  
  list(
    cooks       = d,
    cutoff      = cutoff,
    n_flagged   = sum(d > cutoff),
    flagged_ix  = which(d > cutoff),
    pass        = all(d <= cutoff)
  )
}
#endregion

# dropped function because many points were flagged due to the large dataset, and this is apparently common

# ────────────────────────────────────────────────────────────────────────────────
# 16d. Hosmer–Lemeshow for goodness-of-fit----

#region HOSMER-LEMESHOW CHECK
hoslem_check <- function(fit, g=10) { # g is the number of groups
  # Extract observed outcome robustly; fit$y is often NULL unless glm(..., y=TRUE)
  mf <- tryCatch(stats::model.frame(fit), error = function(e) NULL)
  y_obs <- if (!is.null(mf)) stats::model.response(mf) else fit$y
  if (is.null(y_obs)) stop("Cannot access model response for HL test; refit with y=TRUE or provide data.")
  # Ensure numeric 0/1 for HL test
  y_num <- .as_binary01(y_obs)
  # Use fitted probabilities from the model
  p_hat <- fit$fitted.values
  res <- ResourceSelection::hoslem.test(y_num, p_hat, g=g) # performs the Hosmer-Lemeshow test
  list(
    chisq=res$statistic,
    df=res$parameter,
    p.value=res$p.value,
    pass=(res$p.value>0.05) # null hypothesis is that model fits well; p>0.05 means we do not reject the null
  )
}
#endregion


# ────────────────────────────────────────────────────────────────────────────────
# 16e. Residuals vs. index for independence ----

#region RESIDUALS VS INDEX
plot_resid_index <- function(fit) {
  res <- residuals(fit, type = "deviance") # extract deviance residuals from the fitted model
  df  <- tibble( # initializes a tibble object named df with the specified columns
    index = seq_along(res), # generates a sequence of integers from 1 to the length of the vector res 
    resid = res # assigns the deviance residuals to the column resid
  )

  # plot residuals vs index
  ggplot(df, aes(x = index, y = resid)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_smooth(method = "loess", se = FALSE) +
    labs(
      x     = "Observation index",
      y     = "Deviance residual"
    ) +
    plot_theme_readable()
}
#endregion


# ────────────────────────────────────────────────────────────────────────────────
# 16f. Master assumption wrapper ----

#region CHECK ALL ASSUMPTIONS
check_all_assumptions <- function(fit,
                                  numeric_vars,
                                  vif_thresh = 5,
                                  vif_type   = c("terms", "predictor"),
                                  hl_groups  = 10,
                                  cook_cutoff = NULL ) {
  
  vif_type <- match.arg(vif_type)
  
  # 1) component + residual plots
  cr_plots <- map(numeric_vars, ~ plot_cr_gg(fit, .x))
  names(cr_plots) <- numeric_vars
  
  # 2) VIF
  vif_out <- vif_check(fit, thresh = vif_thresh, type = vif_type)
  
  # 3) Cook's distance (with optional custom cutoff)
  cook_out <- cooks_check(fit, cutoff = cook_cutoff)
  # if no input from user, the cutoff is set to 4/n, where n is the number of observations
  
  # 4) Hosmer–Lemeshow
  hl_out <- hoslem_check(fit, g = hl_groups)
  
  # 5) ROC curve + AUC
  roc_plot <- plot_roc_gg(fit)
  
  # 6) independence: residuals vs. index
  indep_plot <- plot_resid_index(fit)
  
  invisible(list(
    linearity    = cr_plots,
    multicol     = vif_out,
    influencers  = cook_out,
    hosmer       = hl_out,
    roc_auc      = roc_plot,
    independence = indep_plot,
    fit = fit
  ))
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 17. Global plot defaults and high-quality saver ----

# set a readable global ggplot theme and provide a single HQ saver to standardize dpi/size

#region PLOT DEFAULTS AND HQ SAVE
set_plot_defaults <- function(base_size = 14) {
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::theme_set(ggplot2::theme_minimal(base_size = base_size))
    ggplot2::theme_update(
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
  }
  invisible(NULL)
}

save_plot_hq <- function(filename, plot, width = 8, height = 6, dpi = 300, scale = 1) {
  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    units    = "in", # inches
    scale    = scale,
    limitsize = FALSE
  )
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 18. Save component + residual plots ----

# this function saves component+residual plots to PNG files
# automated saving with unique filenames for assumption checking

#region SAVE CR PLOTS
save_crs <- function(checks, suffix, ass_dir) {
    # checks: linearity plots for each numeric variable
    # suffix: file suffix to append (e.g., "_resid")

  purrr::walk2( # iterate over two parallel inputs: names of the elements in checks$linearity and elements in checks$linearity themselves
    names(checks$linearity),
    checks$linearity,
    ~ save_plot_hq( # saves component + residual plots
      filename = file.path(ass_dir, paste0("cr_", .x, suffix, ".png")), # unique filename for each plot
      plot     = .y,
      width    = 8, height = 6, dpi = 300
    )
  )
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 19. Export consolidated results ----

# this function exports streamlined results table with key findings and diagnostics
# combines PRS results with model diagnostics in one interpretable table

#region EXPORT CONSOLIDATED RESULTS
export_consolidated_results <- function(model_list, check_list, output_dir, run_label,
                                       prefix = "results_", alpha = 0.05) {
  
  # Create consolidated results combining PRS effects with diagnostics
  results_df <- purrr::imap_dfr(model_list, function(mod, model_name) {
    
    # Get model results for PRS terms only
    model_results <- broom::tidy(mod) %>%
      filter(grepl("prs_", term, ignore.case = TRUE)) %>%
      mutate(
        model = model_name,
        OR = exp(estimate),
        LCL = exp(estimate - 1.96 * std.error),
        UCL = exp(estimate + 1.96 * std.error),
        p.adj = p.adjust(p.value, method = "bonferroni")
      )
    
    # Get diagnostics
    checks <- check_list[[model_name]]
    
    # Add diagnostic information
    model_results <- model_results %>%
      mutate(
  # Model fit
  AIC = AIC(mod),
  Pseudo_R2 = tryCatch(as.numeric(pscl::pR2(mod)["r2CU"]), error = function(e) NA_real_),
        
        # Assumption checks
        VIF_pass = checks$multicol$pass,
        Cooks_pass = checks$influencers$pass,
        Hosmer_Lemeshow_p = round(checks$hosmer$p.value, 4),
        HL_pass = checks$hosmer$pass,
        
        # Multiple testing
        Bonferroni_threshold = alpha / length(model_results$p.value),
        significant = p.adj < (alpha / length(model_results$p.value))
      ) %>%
      select(model, term, OR, LCL, UCL, p.value, p.adj, significant,
             AIC, Pseudo_R2, VIF_pass, Cooks_pass,
             Hosmer_Lemeshow_p, HL_pass, Bonferroni_threshold)
    
    return(model_results)
  })
  
  # Write to CSV with comprehensive header
  fname <- file.path(output_dir, paste0(prefix, run_label, ".csv"))
  
  # detailed header in the same csv file explaining the results
  header_lines <- c(
    paste("# PRS-Migraine Association Results -", Sys.Date()),
    "#",
    "# KEY FINDINGS:",
    "# OR = Odds Ratio (effect size)",
    "# LCL/UCL = 95% Confidence Interval bounds", 
    "# p.value = uncorrected p-value",
    "# p.adj = Bonferroni-corrected p-value",
    "# significant = TRUE if p.adj < Bonferroni threshold",
    "#",
    "# MODEL DIAGNOSTICS:",
    "# AIC = Akaike Information Criterion (lower = better fit)",
  "# Pseudo_R2 = Nagelkerke's pseudo R-squared (scaled to [0,1])",
    "# VIF_pass = TRUE if all VIF < 5 (no multicollinearity issues)",
    "# VIF_max = Maximum variance inflation factor",
    "# Cooks_pass = TRUE if no influential outliers detected",
    "# Cooks_flagged = Number of observations flagged by Cook's distance",
    "# Hosmer_Lemeshow_p = Goodness-of-fit test p-value (>0.05 = good fit)",
    "# HL_pass = TRUE if Hosmer-Lemeshow test passed",
    "#",
    "# MODELS TESTED:",
    "# noInt_with/without = No interaction model (with/without outliers)",
    "# prsprs_with/without = PRS×PRS interaction model (with/without outliers)",
    "# sex_with/without = Sex interaction model (with/without outliers)", 
    "# both_with/without = All interactions model (with/without outliers)",
    "#"
  )
  
  # Write header and data
  writeLines(header_lines, fname)
  
  # Convert results to character matrix to append with column names
  write.table(results_df, fname, append = TRUE, sep = ",", 
              row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  invisible(results_df)
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 20. K-fold cross-validation AUC ----
# AI transparency: this function was heavily reviewed, improved and debugged with the help of Copilot

# this function performs k-fold cross-validation with AUC calculation and error handling
# comprehensive model validation with automated error handling

#region CV AUC
cv_auc <- function(df, formula_str, k = 5, seed = 42) { # 42: It’s arbitrary; a common convention (Hitchhiker’s “answer”). It’s not statistically special
  set.seed(seed)

  outcome_var <- all.vars(as.formula(formula_str))[1] # extracts the name of the outcome (1st item) represented as a string in the variable formula_str

  y <- df[[outcome_var]] # extracts the outcome variable values from the dataframe
  
  # ensure y is numeric 0/1 robustly
  if (!is.numeric(y) || length(unique(na.omit(y))) > 2 || !all(na.omit(y) %in% c(0,1))) { # if the values are not numeric and the unique values are not 
                                                                                          # more than 2 and the values are not all 0 or 1
    y <- .as_binary01(y) # convert to binary 0/1 using previously defined function
  }

  unique_vals <- sort(unique(y[!is.na(y)])) # gets unique values of y
  cat("Debug: Unique outcome values:", unique_vals, "\n") # prints unique outcome values

  if (!all(unique_vals %in% c(0, 1))) { # if unique values are not 0 or 1, stops and throws an error
    stop("Outcome variable must have values that can be converted to 0 and 1. Found: ",
         paste(unique_vals, collapse = ", ")) # pasting the actual unique values
  }

  y_fac <- factor(y, levels = c(0,1), labels = c("control","case")) #  transforms a outcome into a factor with specified levels and labels

  case_prop <- mean(y == 1, na.rm = TRUE) # calculates the proportion of cases
  cat("Debug: Case proportion:", round(case_prop, 4), "\n") # prints case proportion

  if (case_prop < 0.01 || case_prop > 0.99) { # if case proportion is very imbalanced raises a warning
    warning("Very imbalanced dataset (case proportion: ", round(case_prop, 4), 
            "). Consider using stratified sampling or adjusting k.")
  }

  complete_rows <- !is.na(y) # identifies complete cases where outcome is not NA
  if (!all(complete_rows)) { # if there are any incomplete cases, removes them and prints a message
    cat("Removing", sum(!complete_rows), "rows with missing outcome\n")
    df <- df[complete_rows, ] # subsets the data frame df, keeping only the rows where complete_rows is TRUE
    y <- y[complete_rows] # subsets binary outcome based on the logical vector complete_rows, keeping only the elements where complete_rows is TRUE
    y_fac <- y_fac[complete_rows] # subsets the factorised outcome in the same way, keeping only the elements where complete_rows is TRUE
  }
  
  # use stratified folds to ensure each fold has both classes
  folds <- caret::createFolds(y_fac, k = k, returnTrain = FALSE)
  
  # Debug: Check fold composition
  for (i in seq_along(folds)) {
    fold_y <- y_fac[folds[[i]]]
    fold_table <- table(fold_y)
    cat("Debug: Fold", i, "composition:", paste(names(fold_table), fold_table, collapse = ", "), "\n")
  }
  
  p_hat <- rep(NA_real_, nrow(df)) # numeric vector of the same length as the number of rows in df, with all entries set to NA
  aucs  <- c() # empty numeric vector

  for (i in seq_along(folds)) { # for each fold in folds
    ix <- folds[[i]] # get the indices for the current fold
    fold_y <- y_fac[ix] # get the factor levels for the current fold
    fold_levels <- unique(fold_y) # get the unique levels for the current fold

    if (length(fold_levels) < 2) { # if < 2 unique levels are present raises a warning
      warning("Fold ", i, " contains only one class (", paste(fold_levels, collapse = ", "), 
              ") - skipping AUC calculation for this fold")
      aucs <- c(aucs, NA_real_) # append NA to AUCs for this fold
      next # skip to the next fold
    }
    
    # fit model on training data (all data except current fold)
    train_data <- df[-ix, ] # subset the data frame df, keeping only the rows where ix is FALSE
    test_data <- df[ix, ] # subset the data frame df, keeping only the rows where ix is TRUE

    tryCatch({ # tryCatch block to handle the case in which the model fitting or prediction fails
      fit <- glm(as.formula(formula_str), data = train_data, family = binomial) # fit the model

      p_fold <- .predict_case_prob(fit, newdata = test_data) # predict P(case) on test fold with function previously defined
      p_hat[ix] <- p_fold # store predictions for the current fold

      rocobj <- pROC::roc(fold_y, p_fold, levels = c("control","case"), direction = "<", quiet = TRUE)
      # creates and stores an ROC curve object in the variable rocobj. This object includes various metrics and components, 
      # such as true positive rates (TPR), false positive rates (FPR), area under the curve (AUC), thresholds, and more 
      aucs <- c(aucs, as.numeric(pROC::auc(rocobj))) # retrieve AUC for the current fold
    }, error = function(e) { # if an error occurs raises a warning
      warning("Error in fold ", i, ": ", e$message)
      aucs <<- c(aucs, NA_real_) # append NA to AUCs for this fold
    })
  }

  list(y = y, y_fac = y_fac, p_hat = p_hat, aucs = aucs) # return a list of results
}
#endregion

# ────────────────────────────────────────────────────────────────────────────────
# 21. Plot ROC from GLM fit ----

# this function creates ROC curve plot directly from GLM fit object  
# builds ROC curve and computes AUC for assumption checking

#region PLOT ROC GG
plot_roc_gg <- function(fit) {  # fit: a glm binomial model
  # ensure to compute predictions on the exact rows used for fitting
  # this avoids misalignment or NA issues that can appear with newdata built differently
  mf <- stats::model.frame(fit) # design frame used in the model
  y  <- stats::model.response(mf) # observed outcome used during fit
  p_case <- .predict_case_prob(fit, newdata = mf) # predicted Pr(case) for those rows

  # normalize outcome to explicit control/case labels expected by pROC::roc
  # - If numeric (0/1), map 0->control, 1->case
  # - If factor, infer which level is the 'case' using our helper, then order as control, case
  if (is.numeric(y)) {
    y_fac <- factor(y, levels = c(0, 1), labels = c("control", "case"))
  } else {
    yf <- factor(y)
    if (nlevels(yf) != 2) stop("Outcome must have exactly 2 levels for ROC plotting")
    lv <- levels(yf)
    case_idx <- .infer_case_level(lv)           # choose the most likely case level
    # Relevel so that control is the other level and case is the inferred one
    control_label <- lv[setdiff(1:2, case_idx)]
    case_label <- lv[case_idx]
    y_fac <- factor(as.character(yf), levels = c(control_label, case_label), labels = c("control","case"))
  }

  # build ROC object; direction "<" means higher scores indicate the positive class (case)
  rocobj <- pROC::roc(y_fac, p_case, levels = c("control","case"), direction = "<", quiet = TRUE)
  
  # extract ROC coordinates for plotting: FPR = 1 - specificity, TPR = sensitivity
  dfroc <- tibble(
    fpr = 1 - rocobj$specificities, # false positive rate
    tpr = rocobj$sensitivities # true positive rate
  )

  # plot ROC curve
  auc_val <- as.numeric(pROC::auc(rocobj))
  ggplot(dfroc, aes(x = fpr, y = tpr)) +
    geom_line(color = "steelblue", size = 1) + # ROC curve
    geom_abline(linetype = "dashed", color = "gray50") + # random-chance baseline
    labs(
      x = "False Positive Rate",
      y = "True Positive Rate",
      title = sprintf("AUC = %.3f", auc_val)
    ) +
    plot_theme_readable() +
    theme(plot.title = element_text(size = 16))   # slightly larger title for readability
}
#endregion