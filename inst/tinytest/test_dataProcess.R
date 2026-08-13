# Test dataProcess with default parameters ------------------------------------
QuantDataDefault = dataProcess(SRMRawData, use_log_file = FALSE)
QuantDataDefaultLinear = dataProcess(DDARawData, use_log_file = FALSE, 
                                     summaryMethod = "linear")

# Test dataProcess with numberOfCores parameter ----------------------
QuantDataParallel = dataProcess(SRMRawData, use_log_file = FALSE,
                                               numberOfCores = 2)
QuantDataParallelLinear = dataProcess(DDARawData, use_log_file = FALSE, 
                                     summaryMethod = "linear", numberOfCores = 2)

expect_equal(nrow(QuantDataDefault$FeatureLevelData),
             nrow(QuantDataParallel$FeatureLevelData))

expect_equal(nrow(QuantDataDefaultLinear$FeatureLevelData),
             nrow(QuantDataParallelLinear$FeatureLevelData))

# SRMRawData is a label-based experiment: heavy ("H") rows must be preserved
# in FeatureLevelData after dataProcess

# Helper for near-equality comparison of data.tables with numeric tolerance
expect_dt_equal <- function(dt1, dt2, cols, tol = 1e-6, info = "") {
    # Check non-numeric columns exactly
    num_cols <- cols[sapply(dt1[, ..cols], is.numeric)]
    key_cols <- cols[!cols %in% num_cols]
    
    # Sort both tables the same way before comparing
    if (length(key_cols) > 0) {
        data.table::setkeyv(dt1, key_cols)
        data.table::setkeyv(dt2, key_cols)
    }
    
    # Check row count
    if (nrow(dt1) != nrow(dt2)) {
        return(tinytest::expect_true(FALSE, info = paste(info, "(row count mismatch)")))
    }
    
    # Non-numeric exact match
    if (length(key_cols) > 0) {
        tinytest::expect_true(
            data.table::fsetequal(dt1[, ..key_cols], dt2[, ..key_cols]),
            info = paste(info, "(non-numeric columns)")
        )
    }
    
    # Numeric approximate match
    for (col in num_cols) {
        tinytest::expect_true(
            all(abs(dt1[[col]] - dt2[[col]]) < tol, na.rm = TRUE) &&
                identical(is.na(dt1[[col]]), is.na(dt2[[col]])),
            info = paste(info, sprintf("(column: %s)", col))
        )
    }
}

quant_data_srm = readRDS(
    system.file("tinytest/processed_data/quant_data_srm.rds", 
                package = "MSstats")
)

dt1 <- data.table::as.data.table(quant_data_srm$ProteinLevelData)
dt2 <- data.table::as.data.table(QuantDataDefault$ProteinLevelData)
cols <- sort(intersect(names(dt1), names(dt2)))
expect_true(
    data.table::fsetequal(dt1[, ..cols], dt2[, ..cols]),
    info = "dataProcess ProteinLevelData for SRMRawData should be identical to previously saved output"
)
dt1 <- data.table::as.data.table(quant_data_srm$FeatureLevelData)
dt2 <- data.table::as.data.table(QuantDataDefault$FeatureLevelData)
cols <- sort(intersect(names(dt1), names(dt2)))
expect_true(
    data.table::fsetequal(dt1[, ..cols], dt2[, ..cols]),
    info = "dataProcess FeatureLevelData for SRMRawData should be identical to previously saved output"
)

expect_true(
    "L" %in% QuantDataDefault$FeatureLevelData$LABEL,
    info = "SRMRawData FeatureLevelData must contain L rows"
)
expect_true(
    nrow(QuantDataDefault$ProteinLevelData) > 0,
    info = "SRMRawData must produce non-empty ProteinLevelData"
)
expect_true(
    "LogIntensities" %in% colnames(QuantDataDefault$ProteinLevelData),
    info = "regression: ProteinLevelData must have LogIntensities column"
)
expect_true(
    "GROUP" %in% colnames(QuantDataDefault$ProteinLevelData),
    info = "regression: ProteinLevelData must have GROUP column"
)


# DDARawData Regression Testing
quant_data_dda = readRDS(
    system.file("tinytest/processed_data/quant_data_dda.rds", 
                package = "MSstats")
)
QuantDataDefaultDDA = dataProcess(DDARawData, use_log_file = FALSE)

dt1 <- data.table::as.data.table(quant_data_dda$ProteinLevelData)
dt2 <- data.table::as.data.table(QuantDataDefaultDDA$ProteinLevelData)
cols <- sort(intersect(names(dt1), names(dt2)))
expect_dt_equal(dt1[, ..cols], dt2[, ..cols], cols,
                info = "dataProcess ProteinLevelData for DDARawData should be identical to previously saved output"
)

dt1 <- data.table::as.data.table(quant_data_dda$FeatureLevelData)
dt2 <- data.table::as.data.table(QuantDataDefaultDDA$FeatureLevelData)
cols <- sort(intersect(names(dt1), names(dt2)))
expect_dt_equal(dt1[, ..cols], dt2[, ..cols], cols,
                info = "dataProcess FeatureLevelData for DDARawData should be identical to previously saved output"
)


# Test dataProcess with technical replicates & fractions ------------------
msstats_input_fractions_techreps = data.table::fread(
    system.file("tinytest/processed_data/input_techreps_fractions.csv",
                package = "MSstats")
)
QuantDataTechRepsFractions = dataProcess(msstats_input_fractions_techreps, 
                                         use_log_file = FALSE)
expect_true(!is.null(QuantDataTechRepsFractions))
expect_true(nrow(QuantDataTechRepsFractions$FeatureLevelData) > 0)

# Test dataProcess with linear summary method and anomaly scores --------------
msstats_input = data.table::data.table(
    ProteinName = c(rep("Q9UFW8", 10), rep("Q96S19", 15)),
    PeptideSequence = c(rep("AEFEEQNVR", 5), rep("TALYVTPLDR", 5),
                        rep("AFPLAEWQPSDVDQR", 5), rep("ASGLLLER", 5),
                        rep("LowAbundancePeptide", 5)),
    PrecursorCharge = rep(2, 25),
    FragmentIon = rep("y3", 25),
    ProductCharge = rep(1, 25),
    IsotopeLabelType = rep("L", 25),
    Condition = rep(c("A", "A", "A", "B", "B"), 5),
    BioReplicate = rep(seq(1:5), 5),
    Run = rep(paste0("Run", seq(1:5)), 5),
    Intensity = c(1000, 1500, 2000, 2500, 3000,
                  1100, 1600, 2100, 2600, 3100,
                  1200, 1700, 2200, 2000, 1700,
                  1300, 1800, 1200, 2800, 1800,
                  100, 200, 300, 400, 500),
    AnomalyScores = c(rep(0.03, 10), c(0.01,0.01,0.01,0.9,0.8), c(0.01,0.01,0.01,0.7,0.5), rep(0.10, 5))
)

# Must run with linear summarization
expect_error(
    dataProcess(msstats_input, summaryMethod="TMP"),
    pattern = "AnomalyScores column detected in your input columns.*set summaryMethod to 'linear'",
    info = "Should error when AnomalyScores is present but summaryMethod is not 'linear'"
)

# Process data with linear summarization
QuantDataLinearAnomaly = dataProcess(msstats_input,
                                     normalization=FALSE, # no normalization
                                     MBimpute=TRUE, # Use imputation
                                     summaryMethod="linear" # Key MSstats+ parameter
)

expect_true("Variance" %in% colnames(QuantDataLinearAnomaly$ProteinLevelData),
            info = "Variance column should be present in ProteinLevelData")
protein_data = QuantDataLinearAnomaly$ProteinLevelData
expect_true(all(protein_data$Variance > 0, na.rm = TRUE),
            info = "All variance values should be positive")
expect_true(all(is.finite(protein_data$Variance), na.rm = TRUE),
            info = "All variance values should be finite")
variance_values = protein_data[protein_data$Protein == "Q9UFW8"]$Variance
expect_true(all(variance_values == variance_values[1]),
            info = paste("Not all variance values are equal for Q9UFW8"))

# Create comparison dataset with low anomaly scores for contrast
msstats_input_low_anomaly = data.table::data.table(
    ProteinName = c(rep("Q9UFW8", 10), rep("Q96S19", 15)),
    PeptideSequence = c(rep("AEFEEQNVR", 5), rep("TALYVTPLDR", 5),
                        rep("AFPLAEWQPSDVDQR", 5), rep("ASGLLLER", 5),
                        rep("LowAbundancePeptide", 5)),
    PrecursorCharge = rep(2, 25),
    FragmentIon = rep("y3", 25),
    ProductCharge = rep(1, 25),
    IsotopeLabelType = rep("L", 25),
    Condition = rep(c("A", "A", "A", "B", "B"), 5),
    BioReplicate = rep(seq(1:5), 5),
    Run = rep(paste0("Run", seq(1:5)), 5),
    Intensity = c(1000, 1500, 2000, 2500, 3000,
                  1100, 1600, 2100, 2600, 3100,
                  1200, 1700, 2200, 2000, 1700,
                  1300, 1800, 1200, 2800, 1800,
                  100, 200, 300, 400, 500),
    AnomalyScores = rep(0.01, 25) # All low anomaly scores
)

QuantDataLowAnomaly = dataProcess(msstats_input_low_anomaly,
                                  normalization=FALSE,
                                  MBimpute=TRUE,
                                  summaryMethod="linear")
high_anomaly_variance = protein_data$Variance
low_anomaly_variance = QuantDataLowAnomaly$ProteinLevelData$Variance
expect_true(mean(high_anomaly_variance, na.rm = TRUE) > mean(low_anomaly_variance, na.rm = TRUE),
            info = "Mean variance should be higher for dataset with high anomaly scores")


# Verify that without AnomalyScores, different behavior occurs
msstats_input_no_anomaly = msstats_input[, !c("AnomalyScores")]
QuantDataNoAnomaly = dataProcess(msstats_input_no_anomaly,
                                 normalization=FALSE,
                                 MBimpute=TRUE,
                                 summaryMethod="linear")
no_anomaly_variance = QuantDataNoAnomaly$ProteinLevelData$Variance
expect_false(identical(high_anomaly_variance, no_anomaly_variance),
             info = "Variance should be different when AnomalyScores are present vs absent")



# Test MSstatsSummarizeSingleLinear function ------------------------------

single_protein = data.table::data.table(
    PROTEIN = c( "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19" ),
    PEPTIDE = c( "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2" ),
    TRANSITION = c( "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1" ),
    FEATURE = c( "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1" ),
    LABEL = c( "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L" ),
    GROUP_ORIGINAL = c( "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B" ),
    SUBJECT_ORIGINAL = c( 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5 ),
    RUN = c( "1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5" ),
    GROUP = c( "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B" ),
    SUBJECT = c( "1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5" ),
    FRACTION = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ),
    INTENSITY = c( 1200, 1300, 100, 1700, 1800, 200, 2200, 1200, 300, 2000, 2800, 400, 1700, 1800, 500 ),
    ANOMALYSCORES = c( 0.01, 0.01, 0.1, 0.01, 0.01, 0.1, 0.01, 0.01, 0.1, 0.9, 0.3, 0.1, 0.8, 0.2, 0.1 ),
    ABUNDANCE = c( 10.229, 10.344, 6.644, 10.731, 10.814, 7.644, 11.103, 10.229, 8.229, 10.966, 11.451, 8.644, 10.731, 10.814, 8.966 ),
    originalRUN = c( "Run1", "Run1", "Run1", "Run2", "Run2", "Run2", "Run3", "Run3", "Run3", "Run4", "Run4", "Run4", "Run5", "Run5", "Run5" ),
    censored = c( FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE ),
    newABUNDANCE = c( 10.229, 10.344, NA, 10.731, 10.814, NA, 11.103, 10.229, NA, 10.966, 11.451, NA, 10.731, 10.814, NA ),
    cen = c( 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0 ),
    nonmissing = c( TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE ),
    n_obs = c( 5, 5, 0, 5, 5, 0, 5, 5, 0, 5, 5, 0, 5, 5, 0 ),
    n_obs_run = c( 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ),
    total_features = c( 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ),
    prop_features = c( 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667 ),
    nonmissing_all = c( TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE ),
    ABUNDANCE_cut = c( 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA ),
    any_censored = c( FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE )
)

summarized = MSstats::MSstatsSummarizeSingleLinear(
    single_protein,
    impute = FALSE,
    censored_symbol = "NA", 
    remove50missing = FALSE, 
    aft_iterations = 90,
    equal_variances = TRUE
)

protein_level_summary = summarized[[1]]
feature_level_summary = summarized[[2]]

## Basic tests
expect_equal(nrow(protein_level_summary), 5)
expect_equal(ncol(protein_level_summary), 5)  # Protein, RUN, LogIntensities, Variance, LABEL
expect_true("LABEL" %in% colnames(protein_level_summary),
            info = "MSstatsSummarizeSingleLinear: result must include LABEL column")
expect_true(all(protein_level_summary$Protein == "Q96S19"))
expect_equal(as.numeric(protein_level_summary$RUN), 1:5)
expect_equal(nrow(feature_level_summary), 10)
expect_equal(ncol(feature_level_summary), 6)  # newABUNDANCE, cen, RUN, FEATURE, LABEL, predicted
expect_equal(length(unique(feature_level_summary$FEATURE)), 2)
expect_equal(length(unique(feature_level_summary$RUN)), 5)

## Verify that variance is higher in runs 4-5 with high anomaly scores
variance_runs_1_3 = protein_level_summary[RUN %in% 1:3, Variance]
variance_runs_4_5 = protein_level_summary[RUN %in% 4:5, Variance]

for (var_45 in variance_runs_4_5) {
    for (var_13 in variance_runs_1_3) {
        expect_true(var_45 > var_13, 
                  info = paste("Variance in runs 4-5 should be > variance in runs 1-3:",
                               var_45, ">", var_13))
    }
}

mean_variance_1_3 = mean(variance_runs_1_3)
mean_variance_4_5 = mean(variance_runs_4_5)

expect_true(mean_variance_4_5 > 10 * mean_variance_1_3,
          info = paste("Mean variance of runs 4-5 should be at least 10x higher:",
                       mean_variance_4_5, "vs", 10 * mean_variance_1_3))

## Verify that protein levels in runs 4-5 are closer to ASGLLLER than AFPLAEWQPSDVDQR
protein_runs_4_5 = protein_level_summary[RUN %in% 4:5]
afpla_runs_4_5 = feature_level_summary[FEATURE == "AFPLAEWQPSDVDQR_2_y3_1" & RUN %in% 4:5]
asgll_runs_4_5 = feature_level_summary[FEATURE == "ASGLLLER_2_y3_1" & RUN %in% 4:5]
for (i in 1:nrow(protein_runs_4_5)) {
    run_num = protein_runs_4_5$RUN[i]
    protein_logint = protein_runs_4_5$LogIntensities[i]
    
    afpla_abundance = afpla_runs_4_5[RUN == run_num, newABUNDANCE]
    asgll_abundance = asgll_runs_4_5[RUN == run_num, newABUNDANCE]
    
    dist_to_afpla = abs(protein_logint - afpla_abundance)
    dist_to_asgll = abs(protein_logint - asgll_abundance)
    
    expect_true(dist_to_asgll < dist_to_afpla,
              info = paste("Run", run_num, ": Protein LogIntensity", protein_logint,
                           "should be closer to ASGLLLER", asgll_abundance,
                           "(dist =", round(dist_to_asgll, 4), ")",
                           "than to AFPLAEWQPSDVDQR", afpla_abundance,
                           "(dist =", round(dist_to_afpla, 4), ")"))
}


# MSstatsSummarizeSingleTMP tests ------------------------------------------

MSstatsConvert::MSstatsLogsSettings(FALSE)

single_protein_labeled <- data.table::data.table(
    PROTEIN       = rep("P1", 4),
    FEATURE       = rep("F1", 4),
    LABEL         = c("H",   "L",   "H",   "L"),
    RUN           = c("R1",  "R1",  "R2",  "R2"),
    newABUNDANCE  = c(10.0,  14.0,  12.0,  16.0),
    is_labeled_ref = c(TRUE,  FALSE, TRUE,  FALSE),
    censored      = rep(FALSE, 4),
    n_obs         = rep(2L,  4),
    n_obs_run     = rep(1L,  4),
    prop_features = rep(1.0, 4)
)

result_labeled <- MSstatsSummarizeSingleTMP(
    single_protein_labeled,
    impute          = FALSE,
    censored_symbol = NULL,
    remove50missing = FALSE,
    aft_iterations  = 90
)

protein_level_labeled <- result_labeled[[1]]

expect_equal(
    sort(as.numeric(protein_level_labeled$LogIntensities)),
    c(15, 15),
    info = paste(
        "MSstatsSummarizeSingleTMP: 'is_labeled_ref' column present with TRUE values must set",
        "is_labeled_reference=TRUE, causing .runTukey to call .adjustLRuns;",
        "both runs should converge to H-adjusted LogIntensities of 15,",
        "not the raw L values (14, 16) that the unlabeled path would return"
    )
)

# MSstatsSummarizeSingleTMP: SRM imputation — H rows must NOT be imputed ------
# For SRM experiments, H is the normalization reference and must never be
# imputed. Only censored L rows (is_labeled_ref=FALSE) should receive a
# predicted value from the survival model.

make_srm_impute_input <- function() {
    runs   <- paste0("R", 1:4)
    levels_rc <- c("0", runs)
    f1 <- data.table::data.table(
        PROTEIN  = "P1",
        FEATURE  = "F1",
        LABEL    = c("H","H","H","H", "L","L","L","L"),
        RUN      = c(runs, runs),
        # F1-H-R1 censored (H reference — must NOT be imputed)
        # F1-L-R2 censored (light peptide — MUST be imputed)
        newABUNDANCE = c(NA,   10.5, 11.0, 11.5,  14.0, NA,   15.0, 15.5),
        censored     = c(TRUE, FALSE,FALSE,FALSE,  FALSE,TRUE, FALSE,FALSE),
        cen          = c(0L,   1L,   1L,   1L,     1L,   0L,   1L,   1L),
        is_labeled_ref = c(TRUE,TRUE,TRUE,TRUE, FALSE,FALSE,FALSE,FALSE)
    )
    f2 <- data.table::data.table(
        PROTEIN  = "P1",
        FEATURE  = "F2",
        LABEL    = c("H","H","H","H", "L","L","L","L"),
        RUN      = c(runs, runs),
        newABUNDANCE = c(10.0,10.5,11.0,11.5, 14.0,14.5,15.0,15.5),
        censored     = rep(FALSE, 8),
        cen          = rep(1L, 8),
        is_labeled_ref = c(TRUE,TRUE,TRUE,TRUE, FALSE,FALSE,FALSE,FALSE)
    )
    dt <- data.table::rbindlist(list(f1, f2))
    dt[, ref_covariate := factor(
        ifelse(is_labeled_ref == FALSE, as.character(RUN), "0"),
        levels = levels_rc
    )]
    dt[, FEATURE   := factor(FEATURE)]
    dt[, RUN       := factor(RUN)]
    dt[, n_obs     := 4L]
    dt[, n_obs_run := 2L]
    dt[, ANOMALYSCORES := NA_real_]
    dt
}

result_srm_imp <- MSstatsSummarizeSingleTMP(
    make_srm_impute_input(),
    impute          = TRUE,
    censored_symbol = "NA",
    remove50missing = FALSE,
    aft_iterations  = 90
)

survival_srm <- result_srm_imp[[2]]

# Censored H reference row: predicted must remain NA (not imputed)
h_cens_pred <- survival_srm[
    as.character(FEATURE) == "F1" &
    as.character(LABEL)   == "H" &
    as.character(RUN)     == "R1",
    predicted
]
expect_true(
    length(h_cens_pred) > 0 && all(is.na(h_cens_pred)),
    info = "MSstatsSummarizeSingleTMP SRM: censored H rows must NOT receive an imputed predicted value"
)

# Censored L row: predicted must be a finite imputed value
l_cens_pred <- survival_srm[
    as.character(FEATURE) == "F1" &
    as.character(LABEL)   == "L" &
    as.character(RUN)     == "R2",
    predicted
]
expect_true(
    length(l_cens_pred) > 0 && all(is.finite(l_cens_pred)),
    info = "MSstatsSummarizeSingleTMP SRM: censored L rows must receive a finite imputed predicted value"
)

# --- Same SRM imputation, but via aft_solver = "cg" ------------------------
# Same invariants must hold (H never imputed, L gets a finite prediction),
# and the imputed values themselves should closely match the default
# aft_solver = "cholesky" path, since both solve the same Newton step.
#
# make_srm_impute_input()'s uncensored values are an exactly noise-free
# linear function of RUN, which makes the Gaussian scale MLE degenerate
# (unbounded as residuals -> 0). That's fine for the qualitative H/L
# invariant checks above, but not a meaningful numeric comparison between
# solvers, so a little jitter is added here to make the fit well-posed.

make_srm_impute_input_with_noise <- function(seed) {
    input <- make_srm_impute_input()
    set.seed(seed)
    input[cen == 1L,
          newABUNDANCE := newABUNDANCE + rnorm(.N, sd = 0.01)]
    input
}

result_srm_imp_chol_noisy <- MSstatsSummarizeSingleTMP(
    make_srm_impute_input_with_noise(seed = 1),
    impute          = TRUE,
    censored_symbol = "NA",
    remove50missing = FALSE,
    aft_iterations  = 90,
    aft_solver      = "cholesky"
)
result_srm_imp_cg <- MSstatsSummarizeSingleTMP(
    make_srm_impute_input_with_noise(seed = 1),
    impute          = TRUE,
    censored_symbol = "NA",
    remove50missing = FALSE,
    aft_iterations  = 90,
    aft_solver      = "cg"
)

survival_srm_chol_noisy <- result_srm_imp_chol_noisy[[2]]
survival_srm_cg <- result_srm_imp_cg[[2]]

h_cens_pred_cg <- survival_srm_cg[
    as.character(FEATURE) == "F1" &
    as.character(LABEL)   == "H" &
    as.character(RUN)     == "R1",
    predicted
]
expect_true(
    length(h_cens_pred_cg) > 0 && all(is.na(h_cens_pred_cg)),
    info = "MSstatsSummarizeSingleTMP SRM (aft_solver = cg): censored H rows must NOT receive an imputed predicted value"
)

l_cens_pred_cg <- survival_srm_cg[
    as.character(FEATURE) == "F1" &
    as.character(LABEL)   == "L" &
    as.character(RUN)     == "R2",
    predicted
]
l_cens_pred_chol_noisy <- survival_srm_chol_noisy[
    as.character(FEATURE) == "F1" &
    as.character(LABEL)   == "L" &
    as.character(RUN)     == "R2",
    predicted
]
expect_true(
    length(l_cens_pred_cg) > 0 && all(is.finite(l_cens_pred_cg)),
    info = "MSstatsSummarizeSingleTMP SRM (aft_solver = cg): censored L rows must receive a finite imputed predicted value"
)
expect_equal(
    l_cens_pred_cg, l_cens_pred_chol_noisy, tolerance = 1e-4,
    check.attributes = FALSE,
    info = "MSstatsSummarizeSingleTMP SRM: aft_solver = cg should closely match aft_solver = cholesky"
)
