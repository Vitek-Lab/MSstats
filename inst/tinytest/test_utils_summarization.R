single_protein <- data.table::data.table(
    PROTEIN = c( "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19" ),
    PEPTIDE = c( "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2" ),
    TRANSITION = c( "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1" ),
    FEATURE = c( "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1" ),
    LABEL = c( "L", "L", "L", "L", "L", "L", "L", "L", "L", "L" ),
    GROUP_ORIGINAL = c( "A", "A", "A", "A", "A", "A", "B", "B", "B", "B" ),
    SUBJECT_ORIGINAL = c( 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 ),
    RUN = c( 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 ),
    GROUP = c( "A", "A", "A", "A", "A", "A", "B", "B", "B", "B" ),
    SUBJECT = c( 1, 1, 2, 2, 3, 3, 4, 4, 5, 5 ),
    FRACTION = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ),
    INTENSITY = c( 1200, 1300, 1700, 1800, 2200, 1200, 2000, 2800, 1700, 1800 ),
    ANOMALYSCORES = c( 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.9, 0.7, 0.8, 0.5 ),
    ABUNDANCE = c( 10.229, 10.344, 10.731, 10.814, 11.103, 10.229, 10.966, 11.451, 10.731, 10.814 ),
    originalRUN = c( "Run1", "Run1", "Run2", "Run2", "Run3", "Run3", "Run4", "Run4", "Run5", "Run5" ),
    censored = c( FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE ),
    newABUNDANCE = c( 10.229, 10.344, 10.731, 10.814, 11.103, 10.229, 10.966, 11.451, 10.731, 10.814 ),
    cen = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ),
    nonmissing = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE ),
    n_obs = c( 5, 5, 5, 5, 5, 5, 5, 5, 5, 5 ),
    n_obs_run = c( 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ),
    total_features = c( 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ),
    prop_features = c( 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667 ),
    nonmissing_all = c( TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE ),
    ABUNDANCE_cut = c( 10.127, 10.127, 10.127, 10.127, 10.127, 10.127, 10.127, 10.127, 10.127, 10.127 ),
    any_censored = c( FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE ),
    weights = c( 100, 100, 100, 100, 100, 100, 1.111, 1.429, 1.25, 2 )
)

summarized = MSstats:::.fitLinearModel(single_protein, FALSE, FALSE, TRUE)

expect_equal(summarized$weights, single_protein$weights,
             info = "Weights from model should match input weights")

single_protein$weights = NULL
summarized_no_weights = MSstats:::.fitLinearModel(single_protein, FALSE, FALSE, TRUE)
expect_true(is.null(summarized_no_weights$weights),
            info = "Weights should be null when no weights are provided")

feature_coefficient_with_weights = 
    summarized$coefficients[['FEATUREASGLLLER_2_y3_1']]
feature_coefficient_no_weights = 
    summarized_no_weights$coefficients[['FEATUREASGLLLER_2_y3_1']]
expect_false(feature_coefficient_with_weights == feature_coefficient_no_weights,
             info = "Feature coefficients should differ when weights are used")


make_srm_input <- function() {
    dt <- data.table::data.table(
        PROTEIN       = rep("P1", 12),
        FEATURE       = factor(rep(c("F1", "F2", "F3"), each = 4)),
        RUN           = factor(rep(c("R1","R2","R3","R4"), 3)),
        LABEL         = rep(c("H","H","L","L"), 3),
        ABUNDANCE     = c(10, 11, 14, 15,
                          10.2, 11.2, 14.2, 15.2,
                          9.8, 10.8, 13.8, 14.8),
        newABUNDANCE  = c(10, 11, 14, 15,
                          10.2, 11.2, 14.2, 15.2,
                          9.8, 10.8, 13.8, 14.8),
        weights       = rep(NA_real_, 12)
    )
    dt[["ref_covariate"]] <- factor(
        rep(c("0","0","R3","R4"), 3),
        levels = c("0","R3","R4")
    )
    dt
}
lm_fit <- MSstats:::.fitLinearModel(
    make_srm_input(),
    is_single_feature = FALSE,
    is_labeled        = TRUE,
    equal_variances   = TRUE
)
expect_true(
    any(grepl("ref_covariate", names(coef(lm_fit)))),
    info = "ref_covariate: model coefficients must include 'ref_covariate' terms"
)
