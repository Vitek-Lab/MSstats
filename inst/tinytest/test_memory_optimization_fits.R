# Tests for memory optimizations in the dataProcess pipeline.
# These tests verify that model objects are stripped of heavy fields (Issue 1).


# --- Issue 1: Model Objects Retain Entire Datasets ----------------------------
#
# .fitSurvival() and .fitLinearModel() used to return full model objects that
# stored copies of the input data, QR decompositions, residual vectors, etc.
# We now strip the heavy fields before returning. These tests verify the
# stripped fields are NULL and that the objects still work for their intended
# downstream use (predict, summary, vcov).

# Test 1a: .fitSurvival strips $y and $linear.predictors ---
#
# survreg() stores $y (the Surv response object) and $linear.predictors
# (a vector of length nrow). Neither is needed by predict.survreg().
# After stripping, object.size should be measurably smaller.

surv_input = data.table::data.table(
    newABUNDANCE = c(10.1, 11.2, 9.5, 10.8, 12.0, 9.0, 11.5, 10.3,
                     10.5, 11.0, 9.8, 10.2, 12.2, 9.3, 11.8, 10.6),
    cen = c(1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1),
    RUN = factor(rep(1:4, each = 4)),
    FEATURE = factor(rep(c("feat1", "feat2", "feat3", "feat4"), 4))
)
# cen = 1 means observed (real measurement), cen = 0 means left-censored
# (below the instrument's detection limit — true value is unknown but at
# most this threshold). survreg uses Surv(..., type='left') to handle this:
# observed rows contribute their exact value to the likelihood, censored
# rows contribute P(true_value <= threshold). We set censored rows to 7.0
# as the upper bound — survreg knows the real value could be lower.
surv_input[cen == 0, newABUNDANCE := 7.0]

surv_fit = MSstats:::.fitSurvival(surv_input, aft_iterations = 90)

# Stripped fields should be NULL
expect_true(is.null(surv_fit$y),
            info = "survreg $y should be stripped to save memory")
expect_true(is.null(surv_fit$linear.predictors),
            info = "survreg $linear.predictors should be stripped to save memory")

# Fields needed by predict() should still exist
expect_false(is.null(surv_fit$coefficients),
             info = "survreg $coefficients must survive stripping")
expect_false(is.null(surv_fit$scale),
             info = "survreg $scale must survive stripping")

# predict() should still work on the stripped object — this is the actual
# downstream use case in MSstatsSummarizeSingleTMP
predictions = predict(surv_fit, newdata = surv_input)
expect_equal(length(predictions), nrow(surv_input),
             info = "predict() must work on the stripped survreg object")
expect_true(all(is.finite(predictions)),
            info = "predict() should return finite values")

# The stripped object should be smaller than an unstripped equivalent.
# Fit the same model without stripping to compare sizes.
unstripped_fit = survival::survreg(
    survival::Surv(newABUNDANCE, cen, type = "left") ~ FEATURE + RUN,
    data = surv_input, dist = "gaussian")
expect_true(object.size(surv_fit) < object.size(unstripped_fit),
            info = paste("Stripped survreg should be smaller.",
                         "Stripped:", object.size(surv_fit),
                         "Unstripped:", object.size(unstripped_fit)))


# Test 1b: .fitLinearModel strips $model ---
#
# lm() stores $model (a full copy of the data used for fitting). This is
# not needed by summary() or vcov(), which only use $qr, $residuals, and
# $coefficients. After stripping, the object should be smaller and
# summary/vcov should still work.

lm_input = data.table::data.table(
    newABUNDANCE = c(10.1, 11.2, 9.5, 10.8, 12.0, 9.0, 11.5, 10.3,
                     10.5, 11.0, 9.8, 10.2, 12.2, 9.3, 11.8, 10.6),
    RUN = factor(rep(1:4, each = 4)),
    FEATURE = factor(rep(c("feat1", "feat2", "feat3", "feat4"), 4)),
    weights = NA
)

lm_fit = MSstats:::.fitLinearModel(lm_input, is_single_feature = FALSE,
                                    is_labeled = FALSE, equal_variances = TRUE)

# $model should be stripped
expect_true(is.null(lm_fit$model),
            info = "lm $model should be stripped to save memory")

# Fields needed by summary() and vcov() should still exist
expect_false(is.null(lm_fit$coefficients),
             info = "lm $coefficients must survive stripping")
expect_false(is.null(lm_fit$qr),
             info = "lm $qr must survive stripping (needed by summary and vcov)")
expect_false(is.null(lm_fit$residuals),
             info = "lm $residuals must survive stripping (needed by summary)")

# summary() and vcov() should still work — these are the actual downstream
# use cases in MSstatsSummarizeSingleLinear
lm_summary = summary(lm_fit)
expect_false(is.null(lm_summary$coefficients),
             info = "summary() must work on the stripped lm object")
lm_vcov = vcov(lm_fit)
expect_true(is.matrix(lm_vcov),
            info = "vcov() must work on the stripped lm object")

# The stripped object should be smaller than an unstripped equivalent.
unstripped_lm = lm(newABUNDANCE ~ FEATURE + RUN, data = lm_input)
expect_true(object.size(lm_fit) < object.size(unstripped_lm),
            info = paste("Stripped lm should be smaller.",
                         "Stripped:", object.size(lm_fit),
                         "Unstripped:", object.size(unstripped_lm)))
