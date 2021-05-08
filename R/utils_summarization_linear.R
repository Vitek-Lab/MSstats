#' Fit a linear model
#' @param input data.table
#' @param is_single_feature logical, if TRUE, data has single feature
#' @param is_labeled logical, if TRUE, data comes from a labeled experiment
#' @param equal_variances logical, if TRUE, equal variances are assumed
#' @return lm or merMod
#' @keywords internal
.fitLinearModel = function(input, is_single_feature, is_labeled,
                           equal_variances) {
    if (!is_labeled) {
        if (is_single_feature) {
            linear_model = lm(ABUNDANCE ~ RUN , data = input)
        } else {
            linear_model = lm(ABUNDANCE ~ FEATURE + RUN , data = input)
        }
    } else {
        if (is_single_feature) {
            linear_model = lm(ABUNDANCE ~ RUN + ref , data = input)
        } else {
            linear_model = lm(ABUNDANCE ~ FEATURE + RUN + ref, data = input)
        }
    }
    if (!equal_variances) {
        linear_model = .updateUnequalVariances(input = input, 
                                               fit = linear_model,
                                               num_iter = 1)
    }
    linear_model
}


#' Adjust model for unequal variances
#' @param input data.table
#' @param fit lm
#' @param num_iter number of iterations
#' @importFrom lme4 lmer
#' @importFrom stats loess
#' @return merMod
#' @keywords internal
.updateUnequalVariances = function(input, fit, num_iter) {
    for (i in 1:num_iter) {
        if (i == 1) {
            if (class(fit) == "lm") {
                abs.resids = data.frame(abs.resids = abs(fit$residuals))
                fitted = data.frame(fitted = fit$fitted.values)
            } else {
                abs.resids = data.frame(abs.resids = abs(resid(fit)))
                fitted = data.frame(fitted = fitted(fit))
            }
            input = data.frame(input, 
                               "abs.resids" = abs.resids, 
                               "fitted" = fitted)
        }
        fit.loess = loess(abs.resids ~ fitted, data = input)
        loess.fitted = data.frame(loess.fitted = fitted(fit.loess))
        input = data.frame(input, "loess.fitted" = loess.fitted)
        
        ## loess fitted valuaes are predicted sd
        input$weight = 1 / (input$loess.fitted ^ 2)
        input = input[, -which(colnames(input) %in% "abs.resids")]
        
        ## re-fit using weight
        if (class(fit) == "lm") {
            wls.fit = lm(formula(fit), data = input, weights = weight)
        } else {
            wls.fit = lmer(formula(fit), data = input, weights = weight)
        }
        
        ## lm or lmer
        if (class(wls.fit) == "lm") {
            abs.resids = data.frame(abs.resids = abs(wls.fit$residuals))
        } else {
            abs.resids = data.frame(abs.resids = abs(resid(wls.fit)))
        }
        input = data.frame(input, "abs.resids" = abs.resids)
        
        input = input[, -which(colnames(input) %in% c("loess.fitted", "weight"))]
    }
    
    return(wls.fit)
}


#' Check if data has less than two features
#' @param input data.table
#' @return logical
#' @keywords internal
.checkSingleFeature = function(input) {
    data.table::uniqueN(input$FEATURE) < 2
}
