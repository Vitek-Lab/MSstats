# Tests that .fitSurvival() and .fitLinearModel() return model objects with
# heavy fields stripped, and that downstream predict/summary/vcov still work.


# --- .fitSurvival: $y and $linear.predictors are stripped ---------------------
#
# Neither field is needed by predict.survreg().

surv_input = data.table::data.table(
    newABUNDANCE = c(10.1, 11.2, 9.5, 10.8, 12.0, 9.0, 11.5, 10.3,
                     10.5, 11.0, 9.8, 10.2, 12.2, 9.3, 11.8, 10.6),
    cen = c(1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1),
    RUN = factor(rep(1:4, each = 4)),
    FEATURE = factor(rep(c("feat1", "feat2", "feat3", "feat4"), 4))
)
# cen = 0 marks left-censored rows; set them to the upper-bound threshold.
surv_input[cen == 0, newABUNDANCE := 7.0]

surv_fit = MSstats:::.fitSurvival(surv_input, aft_iterations = 90)

expect_true(is.null(surv_fit$y),
            info = "survreg $y should be stripped")
expect_true(is.null(surv_fit$linear.predictors),
            info = "survreg $linear.predictors should be stripped")

expect_false(is.null(surv_fit$coefficients),
             info = "survreg $coefficients must survive stripping")
expect_false(is.null(surv_fit$scale),
             info = "survreg $scale must survive stripping")

predictions = predict(surv_fit, newdata = surv_input)
expect_equal(length(predictions), nrow(surv_input),
             info = "predict() must work on the stripped survreg object")
expect_true(all(is.finite(predictions)),
            info = "predict() should return finite values")

unstripped_fit = survival::survreg(
    survival::Surv(newABUNDANCE, cen, type = "left") ~ FEATURE + RUN,
    data = surv_input, dist = "gaussian")
expect_true(object.size(surv_fit) < object.size(unstripped_fit),
            info = paste("Stripped survreg should be smaller.",
                         "Stripped:", object.size(surv_fit),
                         "Unstripped:", object.size(unstripped_fit)))


# --- .fitLinearModel: $model is stripped --------------------------------------
#
# $model is not needed by summary() or vcov().

lm_input = data.table::data.table(
    newABUNDANCE = c(10.1, 11.2, 9.5, 10.8, 12.0, 9.0, 11.5, 10.3,
                     10.5, 11.0, 9.8, 10.2, 12.2, 9.3, 11.8, 10.6),
    RUN = factor(rep(1:4, each = 4)),
    FEATURE = factor(rep(c("feat1", "feat2", "feat3", "feat4"), 4)),
    weights = NA
)

lm_fit = MSstats:::.fitLinearModel(lm_input, is_single_feature = FALSE,
                                    is_labeled = FALSE, equal_variances = TRUE)

expect_true(is.null(lm_fit$model),
            info = "lm $model should be stripped")

expect_false(is.null(lm_fit$coefficients),
             info = "lm $coefficients must survive stripping")
expect_false(is.null(lm_fit$qr),
             info = "lm $qr must survive stripping (needed by summary and vcov)")
expect_false(is.null(lm_fit$residuals),
             info = "lm $residuals must survive stripping (needed by summary)")

lm_summary = summary(lm_fit)
expect_false(is.null(lm_summary$coefficients),
             info = "summary() must work on the stripped lm object")
lm_vcov = vcov(lm_fit)
expect_true(is.matrix(lm_vcov),
            info = "vcov() must work on the stripped lm object")

unstripped_lm = lm(newABUNDANCE ~ FEATURE + RUN, data = lm_input)
expect_true(object.size(lm_fit) < object.size(unstripped_lm),
            info = paste("Stripped lm should be smaller.",
                         "Stripped:", object.size(lm_fit),
                         "Unstripped:", object.size(unstripped_lm)))
