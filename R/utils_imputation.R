#' @importFrom data.table uniqueN
#' @importFrom survival survreg
#' @import survival
.fitSurvival = function(input) {
    missingness_filter = !is.na(input$ABUNDANCE)
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
            fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ RUN + ref,
                                    data = input, dist = "gaussian")
        } else {
            if (countdf) {
                fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ RUN + ref,
                                        data = input, dist = "gaussian")
            } else {
                fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type='left') ~ FEATURE + RUN + ref,
                                        data = input, dist = "gaussian")
            }
        }
    } else {
        if (n_features == 1L) {
            fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type = "left") ~ RUN,
                                    data = input, dist = "gaussian")    
        } else {
            if (countdf) {
                fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type = "left") ~ RUN,
                                        data = input, dist = "gaussian")
            } else {
                fit = survival::survreg(survival::Surv(ABUNDANCE, cen, type = "left") ~ FEATURE + RUN,
                                        data = input, dist = "gaussian")
            }
        }  
    }
    fit
}


.addSurvivalPredictions = function(input) {
    if (!all(input$n_obs == 0) & !all(input$n_obs_run <= 1)) {
        survival_fit = .fitSurvival(input[LABEL == "L", ])
        ifelse(input$censored & input$LABEL == "L",
               predict(survival_fit, newdata = input, 
                       type = "response"), input$ABUNDANCE)
    } else {
        NA
    }
}
