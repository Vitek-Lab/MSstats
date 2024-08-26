#' @importFrom data.table uniqueN
#' @importFrom survival survreg Surv
#' @keywords internal
.fitSurvival = function(input) {
    FEATURE = RUN = NULL
    
    missingness_filter = is.finite(input$newABUNDANCE)
    n_total = nrow(input[missingness_filter, ])
    n_features = data.table::uniqueN(input[missingness_filter, FEATURE])
    n_runs = data.table::uniqueN(input[missingness_filter, RUN])
    is_labeled = data.table::uniqueN(input$LABEL) > 1
    countdf = n_total  < n_features + n_runs - 1

    # TODO: set.seed here?
    set.seed(100)
    if (is_labeled) {
        if (length(unique(input$FEATURE)) == 1) {
            # with single feature, not converge, wrong intercept
            # need to check
            fit = survreg(Surv(newABUNDANCE, cen, type='left') ~ RUN + ref,
                                    data = input, dist = "gaussian",
                          control = list(maxiter=90))
        } else {
            if (countdf) {
                fit = survreg(Surv(newABUNDANCE, cen, type='left') ~ RUN + ref,
                                        data = input, dist = "gaussian",
                              control = list(maxiter=90))
            } else {
                fit = survreg(Surv(newABUNDANCE, cen, type='left') ~ FEATURE + RUN + ref,
                                        data = input, dist = "gaussian",
                              control = list(maxiter=90))
            }
        }
    } else {
        if (n_features == 1L) {
            fit = survreg(Surv(newABUNDANCE, cen, type = "left") ~ RUN,
                                    data = input, dist = "gaussian",
                          control = list(maxiter=90))
        } else {
            if (countdf) {
                fit = survreg(Surv(newABUNDANCE, cen, type = "left") ~ RUN,
                                        data = input, dist = "gaussian",
                              control = list(maxiter=90))
            } else {
                fit = survreg(Surv(newABUNDANCE, cen, type = "left") ~ FEATURE + RUN,
                                        data = input, dist = "gaussian", 
                              control = list(maxiter=90))
            }
        }  
    }
    fit
}


#' Get predicted values from a survival model
#' @param input data.table
#' @return numeric vector of predictions
#' @importFrom stats predict
#' @keywords internal
.addSurvivalPredictions = function(input) {
    LABEL = NULL
    
    survival_fit = .fitSurvival(input[LABEL == "L", ])
    predict(survival_fit, newdata = input)
}


.modelFeatureIntensity = function(input){
    input[, mean_feat_int :=mean(ABUNDANCE, na.rm=TRUE), 
          by=list(PROTEIN, FEATURE)]
    
    upper = unname(quantile(input$mean_feat_int, 0.9, na.rm=TRUE))
    input$REASON = ifelse(input$mean_feat_int >= upper, 
                          "MAR", "MNAR")
    input$mean_feat_int = NULL
    return(input)
}