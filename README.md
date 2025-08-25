# MSc-Bioinformatics-project-scripts
Scripts used to obtain the results in the MSc thesis 

Documentation for functions contained in file functions_3030116.

    load_packages = "Function: load_packages()
------------------------------------------------------------
Checks for required CRAN packages (haven, dplyr, ggplot2, broom, purrr,
car, ResourceSelection, pROC, caret, rms, pscl), installs missing ones, and loads them.
@param None
@return None",

    load_ukbb_data = "Function: load_ukbb_data()
------------------------------------------------------------
Loads a UK Biobank .dta/.csv/.txt file (haven::read_dta or read.csv) with basic validation.
@param filepath  string — path to data file
@return Data frame",

    clean_migraine_status = "Function: clean_migraine_status()
------------------------------------------------------------
Filters to rows where migraine_col ∈ {0,1}, drops NAs.
@param df           data.frame
@param migraine_col string — column name for binary outcome
@return Filtered data.frame",

    select_columns = "Function: select_columns()
------------------------------------------------------------
Subsets a data.frame to the user‐specified columns.
@param df           data.frame
@param cols_to_keep character vector — column names to retain
@return Subset data.frame",

    remove_missing_prs = "Function: remove_missing_prs()
------------------------------------------------------------
Drops rows with NA in any of the specified PRS columns.
@param df      data.frame
@param cols_na character vector — PRS column names
@return Filtered data.frame",

    remove_outliers_z = "Function: remove_outliers_z()
------------------------------------------------------------
Removes rows where any numeric variable has |Z| > threshold.
@param df           data.frame
@param numeric_vars character vector — numeric column names
@param thresh       numeric — Z‐score cutoff
@return Filtered data.frame",

    plot_theme_readable = "Function: plot_theme_readable()
------------------------------------------------------------
Sets a readable ggplot2 minimal theme with larger text.
@param sz numeric — base font size
@return ggplot2 theme",

    plot_sex_distribution = "Function: plot_sex_distribution()
------------------------------------------------------------
Creates stacked proportion bars of sex distribution by status.
@param df         data.frame
@param status_col string — binary outcome column
@param sex_col    string — sex column
@return ggplot object",

    plot_chip_distribution = "Function: plot_chip_distribution()
------------------------------------------------------------
Creates stacked proportion bars of genotyping array by status.
@param df         data.frame
@param status_col string — binary outcome column
@param chip_col   string — genotyping array column
@return ggplot object",

    choose_2sample_test = "Function: choose_2sample_test()
------------------------------------------------------------
Chooses Welch t-test by default; falls back to Wilcoxon if severe outliers.
@param x numeric vector
@param y numeric vector
@return List(test, p.value, fit, n_x, n_y)",

    prelimin_eda = "Function: prelimin_eda()
------------------------------------------------------------
Comprehensive preliminary EDA with group comparisons and effect sizes; writes a summary file.
@param df           data.frame
@param status_var   string — binary outcome column
@param sex_var      string — sex column
@param age_var      string — age column
@param chip_var     string — genotyping array column
@param female_code  value — code used for 'Female' in sex_var
@param outfile      string — output filename (written inside output_dir)
@param output_dir   string — directory to write the summary
@param dataset_label string — label to show in header
@param n_with_outliers integer — N for with-outliers dataset to display alongside
@return List with key test results (invisible)",

    plot_prs_decile = "Function: plot_prs_decile()
------------------------------------------------------------
Decile-prevalence plot for a PRS variable.
@param df         data.frame
@param prs_var    string — PRS column name
@param status_var string — outcome column name
@return ggplot object",

    plot_uni_or = "Function: plot_uni_or()
------------------------------------------------------------
Vertical forest plot of unadjusted ORs (95% CI) for covariates.
@param df         data.frame
@param covars     character vector — covariate names
@param status_var string — outcome column name
@param outfile    string|NULL — optional PNG path
@return ggplot object",

    eda_save_prs_decile = "Function: eda_save_prs_decile()
------------------------------------------------------------
Saves PRS decile-prevalence plot to file.
@param df         data.frame
@param prs_var    string — PRS column name
@param status_var string — outcome column name
@param outdir     string — output directory
@return NULL (saves file)",

    eda_save_uni_or = "Function: eda_save_uni_or()
------------------------------------------------------------
Saves univariate OR forest plot to file.
@param df         data.frame
@param covariates character vector — covariate names
@param status_var string — outcome column name
@param outdir     string — output directory
@return NULL (saves file)",

    run_all_eda = "Function: run_all_eda()
------------------------------------------------------------
Runs complete EDA pipeline (violin/box, decile, univariate OR; sex/chip distributions) and saves outputs.
@param ukbb_clean data.frame — main dataset
@param prs_cols   character vector — PRS column names (length >= 2)
@param status_col string — outcome column name
@param covariates character vector — covariate names
@param output_dir string — output directory
@param facet_vars character vector — up to two vars (sex, chip) for extra plots
@param subdir     string — subdirectory name (default 'EDA')
@param viz_z_thresh numeric — visualization-only Z threshold (default 5)
@return Invisible list of output file paths",

    fit_prs_model = "Function: fit_prs_model()
------------------------------------------------------------
Fits logistic regression with two PRS and optional PRS×PRS and PRS×sex interactions; computes ORs and CIs.
@param df             data.frame
@param outcome_var    string
@param prs_vars       character vector (length 2)
@param covariates     character vector
@param interaction    logical — include PRS×PRS interaction
@param sex_var        string — sex variable name
@param sex_interaction logical — include PRS×sex interactions
@param alpha          numeric — significance level
@return List(model, results, pvalues, p.adj, threshold, sig_terms)",

    run_sex_stratified_no_interactions = "Function: run_sex_stratified_no_interactions()
------------------------------------------------------------
Fits separate logistic models by sex without any interactions; returns tidy results per stratum.
@param df             data.frame
@param outcome_var    string
@param prs_vars       character vector
@param covariates     character vector (sex excluded)
@param sex_var        string — sex variable name
@return List(models = list(female = ..., male = ...), results = data.frame)",

    plot_cr_gg = "Function: plot_cr_gg()
------------------------------------------------------------
Component+residual plot for a numeric predictor from a glm fit.
@param fit object — glm
@param var string — predictor name
@return ggplot object",

    vif_check = "Function: vif_check()
------------------------------------------------------------
Variance Inflation Factors for multicollinearity (custom for glm; car::vif fallback otherwise).
@param fit    object — glm or lm
@param thresh numeric — VIF cutoff
@param type   'terms' or 'predictor'
@return List(vif_values, pass)",

    cooks_check = "Function: cooks_check()
------------------------------------------------------------
Cook's distance diagnostics with configurable cutoff (default 1 if unspecified).
@param fit    object — glm
@param cutoff numeric — Cook's distance cutoff
@return List(cooks, cutoff, n_flagged, flagged_ix, pass)",

    hoslem_check = "Function: hoslem_check()
------------------------------------------------------------
Hosmer–Lemeshow goodness‐of‐fit test (ResourceSelection::hoslem.test).
@param fit object — glm
@param g   integer — number of groups
@return List(chisq, df, p.value, pass)",

    plot_resid_index = "Function: plot_resid_index()
------------------------------------------------------------
Residuals vs. index plot for independence check.
@param fit object — glm
@return ggplot object",

    check_all_assumptions = "Function: check_all_assumptions()
------------------------------------------------------------
Runs linearity (CR plots), VIF, Cook's D, Hosmer–Lemeshow, ROC, and residuals-index plots.
@param fit          object — glm
@param numeric_vars character vector
@param vif_thresh   numeric
@param vif_type     char — 'terms'/'predictor'
@param hl_groups    integer
@param cook_cutoff  numeric — Cook's distance cutoff
@return List(linearity, multicol, influencers, hosmer, roc_auc, independence, fit)",

    set_plot_defaults = "Function: set_plot_defaults()
------------------------------------------------------------
Sets a global ggplot2 theme_minimal with a chosen base size.
@param base_size numeric — base font size
@return None (sets theme)",

    save_plot_hq = "Function: save_plot_hq()
------------------------------------------------------------
Saves a ggplot with common high-quality defaults.
@param filename string — output path
@param plot     ggplot object
@param width    numeric — inches
@param height   numeric — inches
@param dpi      integer — resolution
@param scale    numeric — scale factor
@return Invisible filename",

    save_crs = "Function: save_crs()
------------------------------------------------------------
Saves component+residual plots (from check_all_assumptions) to Assumptions subfolder.
@param checks list — output from check_all_assumptions
@param suffix string — filename suffix
@param ass_dir string — output directory
@return NULL (saves files)",

    export_consolidated_results = "Function: export_consolidated_results()
------------------------------------------------------------
Exports streamlined results table with key findings and diagnostics (writes CSV to Results).
@param model_list list — named list of fitted models
@param check_list list — named list of assumption checks
@param output_dir string — output directory
@param run_label  string — run identifier
@param prefix     string — filename prefix (default: 'results_')
@param alpha      numeric — significance level (default: 0.05)
@return Data frame (also saves to CSV)",

    diagnose_cv_data = "Function: diagnose_cv_data()
------------------------------------------------------------
Lightweight diagnostic printout for CV outcome variable (class, range, table, NA/infinite counts).
@param df          data.frame
@param outcome_var string — outcome column name
@return None (prints to console)",

    cv_auc = "Function: cv_auc()
------------------------------------------------------------
K-fold cross-validation with AUC calculation (caret stratified folds, pROC AUC).
@param df          data.frame
@param formula_str string — model formula
@param k           integer — number of folds
@param seed        integer — random seed
@return List(y, y_fac, p_hat, aucs)",

    plot_roc_gg = "Function: plot_roc_gg()
------------------------------------------------------------
Builds ROC curve and computes AUC for a glm fit (pROC); returns the ggplot.
@param fit object — glm
@return ggplot object"
