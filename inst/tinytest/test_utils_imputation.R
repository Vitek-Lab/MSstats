make_surv_labeled_single <- function() {
    runs <- paste0("R", 1:3)
    dt <- data.table::rbindlist(list(
        data.table::data.table(
            FEATURE      = factor(rep("F1", 9)),
            RUN          = factor(rep(runs, each = 3)),
            LABEL        = "H",
            newABUNDANCE = seq(10.1, by = 0.1, length.out = 9),
            cen          = 1L
        ),
        data.table::data.table(
            FEATURE      = factor(rep("F1", 9)),
            RUN          = factor(rep(runs, each = 3)),
            LABEL        = "L",
            newABUNDANCE = seq(14.1, by = 0.1, length.out = 9),
            cen          = 1L
        )
    ))
    ref_vals <- ifelse(dt$LABEL == "L", as.character(dt$RUN), "0")
    dt[["ref_covariate"]] <- factor(ref_vals, levels = c("0", runs))
    dt
}

make_surv_labeled_multi_welldetermined <- function() {
    features <- paste0("F", 1:3)
    runs <- paste0("R", 1:4)
    dt <- data.table::CJ(FEATURE = features, RUN = runs, LABEL = c("H", "L"))
    dt[, FEATURE := factor(FEATURE)]
    dt[, RUN := factor(RUN)]
    dt[, newABUNDANCE := 10 + as.integer(FEATURE) + as.integer(RUN) * 0.5 +
           ifelse(LABEL == "L", 4, 0)]
    dt[, cen := 1L]
    ref_vals <- ifelse(dt$LABEL == "L", as.character(dt$RUN), "0")
    dt[["ref_covariate"]] <- factor(ref_vals, levels = c("0", runs))
    dt
}

make_surv_labeled_underdetermined <- function() {
    runs <- c("R1", "R2", "R3")
    dt_h <- data.table::data.table(
        FEATURE      = factor(paste0("F", 1:8)),
        RUN          = factor(rep_len(runs, 8)),
        LABEL        = "H",
        newABUNDANCE = seq(10, by = 0.5, length.out = 8),
        cen          = 1L
    )
    dt_l <- data.table::data.table(
        FEATURE      = factor("F1"),
        RUN          = factor("R1"),
        LABEL        = "L",
        newABUNDANCE = 14,
        cen          = 1L
    )
    dt <- data.table::rbindlist(list(dt_h, dt_l))
    ref_vals <- ifelse(dt$LABEL == "L", as.character(dt$RUN), "0")
    dt[["ref_covariate"]] <- factor(ref_vals, levels = c("0", runs))
    dt
}

make_surv_unlabeled_single <- function() {
    data.table::data.table(
        FEATURE      = factor(rep("F1", 15)),
        RUN          = factor(rep(paste0("R", 1:5), each = 3)),
        LABEL        = "L",
        newABUNDANCE = seq(10, by = 0.5, length.out = 15),
        cen          = 1L
    )
}

make_surv_unlabeled_multi_welldetermined <- function() {
    dt <- data.table::CJ(
        FEATURE = paste0("F", 1:3),
        RUN     = paste0("R", 1:5)
    )
    dt[, FEATURE := factor(FEATURE)]
    dt[, RUN := factor(RUN)]
    dt[, LABEL := "L"]
    dt[, newABUNDANCE := seq(10, by = 0.5, length.out = .N)]
    dt[, cen := 1L]
    dt
}

coef_names <- function(fit) names(coef(fit))
surv_labeled_single <- MSstats:::.fitSurvival(make_surv_labeled_single(), 90)

expect_true(
    any(grepl("ref_covariate", coef_names(surv_labeled_single))),
    info = ".fitSurvival labeled single-feature: ref_covariate must appear in coefficients"
)
expect_false(
    any(grepl("^FEATURE", coef_names(surv_labeled_single))),
    info = ".fitSurvival labeled single-feature: FEATURE must not appear (only one feature)"
)

surv_labeled_multi_wd <- MSstats:::.fitSurvival(make_surv_labeled_multi_welldetermined(), 90)
expect_true(
    any(grepl("ref_covariate", coef_names(surv_labeled_multi_wd))),
    info = ".fitSurvival labeled multi well-determined: ref_covariate must appear in coefficients"
)
expect_true(
    any(grepl("^FEATURE", coef_names(surv_labeled_multi_wd))),
    info = ".fitSurvival labeled multi well-determined: FEATURE must appear in coefficients"
)

surv_labeled_under <- MSstats:::.fitSurvival(make_surv_labeled_underdetermined(), 90)
expect_true(
    any(grepl("ref_covariate", coef_names(surv_labeled_under))),
    info = ".fitSurvival labeled underdetermined: ref_covariate must appear in fallback coefficients"
)
expect_false(
    any(grepl("^FEATURE", coef_names(surv_labeled_under))),
    info = ".fitSurvival labeled underdetermined: FEATURE must not appear in fallback formula"
)

surv_unlabeled_single <- MSstats:::.fitSurvival(make_surv_unlabeled_single(), 90)
expect_false(
    any(grepl("ref_covariate", coef_names(surv_unlabeled_single))),
    info = ".fitSurvival unlabeled single-feature: ref_covariate must not appear"
)
expect_false(
    any(grepl("^FEATURE", coef_names(surv_unlabeled_single))),
    info = ".fitSurvival unlabeled single-feature: FEATURE must not appear (only one feature)"
)

surv_unlabeled_multi_wd <- MSstats:::.fitSurvival(make_surv_unlabeled_multi_welldetermined(), 90)
expect_false(
    any(grepl("ref_covariate", coef_names(surv_unlabeled_multi_wd))),
    info = ".fitSurvival unlabeled multi well-determined: ref_covariate must not appear"
)
expect_true(
    any(grepl("^FEATURE", coef_names(surv_unlabeled_multi_wd))),
    info = ".fitSurvival unlabeled multi well-determined: FEATURE must appear in coefficients"
)
