# Call flag_noninf_data() function as follows: 
# annotatedData <- flag_noninf_data(processedData)

#############################################
## feature selection module
#############################################

#' @importFrom MASS rlm
#' @import tidyr
#' @import dplyr
#' @import purrr
#' @import stringr
#' @import broom
#' @importFrom stats qbinom
#' @importFrom utils data

## Wrapper to fit robust linear model for one protein

#' @keywords internal
#' 
fit_prot_huber <- function(df_prot) {
    fit <- rlm(log2inty ~ run + feature, data = df_prot, scale.est = "Huber")
    
    return(fit)
}

## Calculate feature variance with data from broom::augment
#' @keywords internal
#' 
calc_fvar <- function(augmented_data, s_resid, rm_olr = FALSE, tol = 3) {
    
    select <- dplyr::select
    
    v_resid <- s_resid ^ 2
    if (rm_olr) {
        augmented_data <- augmented_data %>% 
            mutate(is_olr = abs(.resid / s_resid) > tol) %>% 
            filter(!is_olr)
    }
    varfeature <- augmented_data %>% 
        mutate(resid_null = log2inty - mean(log2inty)) %>% 
        group_by(feature) %>%
        summarise(
            nb_run = n(), 
            svar_feature = sum(.resid ^ 2) / (nb_run - 1) / v_resid, 
            svar_ref = sum(resid_null ^ 2) / (nb_run - 1) / v_resid
        ) %>% 
        select(feature, svar_feature, svar_ref)
    
    return(varfeature)
}

## Flag outlier with data from broom::augment
#' @keywords internal
#' 
flag_outlier <- function(augmented_data, s_resid, tol = 3, keep_run = FALSE) {
    
    select <- dplyr::select
    
    outlier <- augmented_data %>% 
        mutate(is_olr = abs(.resid / s_resid) > tol) %>% 
        select(run, feature, is_olr)
    # To keep runs from being completely removed
    if (keep_run) {
        uncovered_run <- outlier %>% 
            group_by(run) %>% 
            filter(all(is_olr)) %>% 
            ungroup() %>% 
            distinct(run)
        # Handle uncovered runs below
        if (nrow(uncovered_run) > 0) {
            outlier <- bind_rows(
                outlier %>% anti_join(uncovered_run), 
                outlier %>% semi_join(uncovered_run) %>% mutate(is_olr = FALSE)
            )
        }
    }
    
    return(outlier)
}

## Flag uninformative features and outliers
## Annotate processed data with two additional columns
##  - feature_quality: ["Uninformative"|"Informative"]
##  - is_outlier: [TRUE|FALSE] 
##
#' @keywords internal
#' 
flag_noninf_data <- function(processedData) {
    
    select <- dplyr::select
    
    # Convert to the working format
    # df_allftr <- processedData %>%
    #     select(
    #         protein = PROTEIN,
    #         peptide = PEPTIDE,
    #         feature = FEATURE, 
    #         run = originalRUN, 
    #         label = LABEL,
    #         log2inty = ABUNDANCE, 
    #         is_censored = censored
    #     ) %>% 
    #     mutate(
    #         is_obs = !(is.na(log2inty) | is_censored), 
    #         log2inty = ifelse(is_obs, log2inty, NA), 
    #         protein = as.character(protein), 
    #         peptide = as.character(peptide), 
    #         feature = as.character(feature), 
    #         run = as.character(run), 
    #         label = as.character(label)
    #     )
    if ("censored" %in% names(processedData)) {
        df_allftr <- processedData %>%
            select(
                protein = PROTEIN,
                peptide = PEPTIDE,
                feature = FEATURE, 
                run = originalRUN, 
                label = LABEL,
                log2inty = ABUNDANCE, 
                is_censored = censored
            )
    } else {
        df_allftr <- processedData %>%
            select(
                protein = PROTEIN,
                peptide = PEPTIDE,
                feature = FEATURE, 
                run = originalRUN, 
                label = LABEL,
                log2inty = ABUNDANCE
            ) %>% 
            mutate(is_censored = FALSE)
    }
    df_allftr <- df_allftr %>% 
        mutate(
            is_obs = !(is.na(log2inty) | is_censored), 
            log2inty = ifelse(is_obs, log2inty, NA), 
            protein = as.character(protein), 
            peptide = as.character(peptide), 
            feature = as.character(feature), 
            run = as.character(run), 
            label = as.character(label)
        )
    
    # Is the dataset from a labeled experiment?
    is_labeled <- any(df_allftr$label == "H")
    all_label <- unique(df_allftr$label)
    
    # Process H/L channels separately
    df_rmftr <- vector("list", length(all_label))
    df_rmpk <- vector("list", length(all_label))
    for (l in seq_along(all_label)) {
        
        lab <- all_label[l]
        df_onelab <- df_allftr %>% filter(label == lab)
        message(paste0("Analyzing features in ", lab, " channel..."))
        
        # In heavy channel, only positive peaks are used for modeling
        if (lab == "H") {
            df_onelab <- df_onelab %>% filter(log2inty > 0)
        }
        
        # Singleton (unreplicated) or undetectable features
        ftr_unrep <- df_onelab %>% 
            group_by(protein, feature) %>% 
            filter(sum(is_obs) <= 1) %>% 
            ungroup() %>% 
            distinct(protein, feature)
        
        # Remove singleton or undetectable features & uncovered runs
        df_onelab <- df_onelab %>% 
            anti_join(ftr_unrep) %>% 
            group_by(protein, run) %>% 
            filter(any(is_obs)) %>% 
            ungroup()
        
        # Coverage of each feature
        obs_protein <- df_onelab %>% 
            group_by(protein) %>% 
            summarise(
                nb_run = n_distinct(run), 
                nb_feature = n_distinct(feature), 
                nb_obs = sum(is_obs)
            ) %>% 
            mutate(nb_full = nb_run * nb_feature, pi_obs = nb_obs / nb_full) %>% 
            mutate(min_obs = qbinom(0.01, nb_run, pi_obs))
        
        # Nested data frame
        nested_prot <- df_onelab %>% 
            select(protein, run, peptide, feature, log2inty, is_obs) %>% 
            group_by(protein) %>% 
            nest()
        
        nested_prot <- nested_prot %>% 
            left_join(obs_protein %>% select(protein, min_obs))
        
        # Coverage
        message("Identifying low-coverage features...")
        
        list_cover <- vector("list", length = nrow(nested_prot))
        for (i in seq_along(list_cover)) {
            augdata <- nested_prot$data[[i]] %>% 
                group_by(feature) %>% 
                mutate(is_lowcvr = sum(is_obs) < nested_prot$min_obs[i]) %>% 
                ungroup()
            nested_prot$data[[i]] <- augdata
            list_cover[[i]] <- augdata %>% distinct(feature, is_lowcvr)
        }
        nested_prot$cover_feature <- list_cover
        rm(list_cover)
        
        # List of low-coverage features
        ftr_lowcvr <- nested_prot %>% 
            unnest(cover_feature) %>% 
            filter(is_lowcvr) %>% 
            distinct(protein, feature)
        
        # Keep proteins with at least 3 good-coverage features
        is_ftreff <- map_lgl(nested_prot$cover_feature, ~sum(!(.$is_lowcvr)) > 2)
        nested_prot <- nested_prot[is_ftreff, ]
        
        # Robust linear model (NB: non-converged models are not used)
        message("Fitting robust linear models...")
        
        list_sdata <- map(nested_prot$data, ~filter(., !is_lowcvr))
        fit <- list_sdata %>% map(possibly(quietly(fit_prot_huber), NULL))
        is_failed <- map_lgl(fit, is.null)
        nested_prot <- nested_prot[!is_failed, ]
        fit <- fit[!is_failed] %>% transpose()
        is_warned <- !map_lgl(fit$warnings, is_empty)
        nested_prot$rlm_fit <- fit$result
        nested_prot <- nested_prot %>% 
            slice(which(!(is_warned))) %>% 
            mutate(df_resid = map_dbl(rlm_fit, ~summary(.)$df[2])) %>% 
            mutate(s_resid = map_dbl(rlm_fit, ~summary(.)$sigma))
        rm(list = c("list_sdata", "fit", "is_failed", "is_warned"))

        # Extract variances and degrees of freedom
        var_protein <- nested_prot %>% 
            select(protein, s_resid, df_resid) %>% 
            mutate(var_resid = s_resid ^ 2) %>% 
            filter(!is.na(s_resid), df_resid > 0)
        
        # Shrinkage variance estimation with limma
        eb_fit <- limma::squeezeVar(var_protein$var_resid, var_protein$df_resid, robust = TRUE)
        var_protein <- var_protein %>% 
            ungroup() %>%
            mutate(var_resid_eb = eb_fit$var.post) %>% 
            mutate(s_resid_eb = sqrt(var_resid_eb))
        
        nested_prot <- nested_prot %>% 
            left_join(var_protein %>% dplyr::select(protein, s_resid_eb))
        
        # Feature variance and outlier detection
        message("Identifying outliers and calculating feature variances...")
        
        list_var <- vector("list", nrow(nested_prot))
        for (i in seq_along(list_var)) {
            s_eb <- nested_prot$s_resid_eb[i]
            if (!is.na(s_eb) & !all(is.na(nested_prot$rlm_fit[[i]]$residuals))) {
                augdata <- nested_prot$rlm_fit[[i]] %>% augment()
                if (lab == "L") {
                    list_var[[i]] <- calc_fvar(augdata, s_eb, rm_olr = TRUE)
                    nested_prot$data[[i]] <- suppressMessages(
                        nested_prot$data[[i]] %>% 
                            left_join(flag_outlier(augdata, s_eb, keep_run = FALSE)) %>% 
                            mutate(is_olr = ifelse(is.na(is_olr), FALSE, is_olr))
                    )
                } else {
                    # No outlier detection in heavy channel; all measurements for calculation
                    list_var[[i]] <- calc_fvar(augdata, s_eb, rm_olr = FALSE)
                    nested_prot$data[[i]] <- nested_prot$data[[i]] %>% mutate(is_olr = FALSE)
                }
            } else {
                nested_prot$data[[i]] <- nested_prot$data[[i]] %>% mutate(is_olr = FALSE)
            }
        }
        nested_prot$var_feature <- list_var
        rm("list_var")
        
        # Noisy features
        allfvar <- nested_prot %>% 
            filter(!is.na(s_resid_eb)) %>% 
            select(protein, var_feature) %>% 
            unnest(var_feature)
        fvar_cut <- quantile(allfvar$svar_ref, 0.05, na.rm = TRUE)
        ftr_hivar <- allfvar %>% filter(svar_feature > fvar_cut) %>% distinct(protein, feature)
        
        # Combine unreplicated, low-coverage, and noisy features
        df_rmftr[[l]] <- bind_rows(ftr_hivar, ftr_lowcvr, ftr_unrep) %>% 
            select(protein, feature) %>% 
            mutate(label = lab)
        
        # Outliers
        df_rmpk[[l]] <- nested_prot %>% 
            unnest(data) %>% 
            filter(is_olr) %>% 
            select(protein, run, feature) %>% 
            mutate(label = lab)
        
        message(paste0("Outlier detection and feature selection in ", lab, " channel are completed"))
    }
    df_rmftr <- bind_rows(df_rmftr)
    df_rmpk <- bind_rows(df_rmpk)
    
    # Non-positive peaks in heavy channel are considered as outliers
    if (is_labeled) {
        df_rmpk <- bind_rows(
            df_rmpk, 
            df_allftr %>% 
                filter(label == "H", log2inty <= 0) %>% 
                select(protein, run, feature, label)
        )
    }
    
    # Annotate feature and peak qualities with 2 additional columns
    df_rmftr <- df_rmftr %>%
        rename(PROTEIN = protein, FEATURE = feature, LABEL = label)
    df_rmpk <- df_rmpk %>%
        rename(PROTEIN = protein, originalRUN = run, FEATURE = feature, LABEL = label)

    annotatedData <- bind_rows(
        suppressWarnings(processedData %>% semi_join(df_rmftr) %>% mutate(feature_quality = "Uninformative")),
        suppressWarnings(processedData %>% anti_join(df_rmftr) %>% mutate(feature_quality = "Informative"))
    )
    annotatedData <- bind_rows(
        suppressWarnings(annotatedData %>% semi_join(df_rmpk) %>% mutate(is_outlier = TRUE)),
        suppressWarnings(annotatedData %>% anti_join(df_rmpk) %>% mutate(is_outlier = FALSE))
    )
    
    return(annotatedData)
}




## Flag uninformative features and outliers
## Annotate processed data with two additional columns
##  - feature_quality: ["Uninformative"|"Informative"]
##  - is_outlier: [TRUE|FALSE] 
##
#' @keywords internal
#' 
flag_noninf_data_nbftr <- function(processedData) {
    
    select <- dplyr::select
    
    ### [TEST] ###
    min_ftr <- 2
    ### [TEST] ###
    
    # Convert to the working format
    # df_allftr <- processedData %>%
    #     select(
    #         protein = PROTEIN,
    #         peptide = PEPTIDE,
    #         feature = FEATURE, 
    #         run = originalRUN, 
    #         label = LABEL,
    #         log2inty = ABUNDANCE, 
    #         is_censored = censored
    #     ) %>% 
    #     mutate(
    #         is_obs = !(is.na(log2inty) | is_censored), 
    #         log2inty = ifelse(is_obs, log2inty, NA), 
    #         protein = as.character(protein), 
    #         peptide = as.character(peptide), 
    #         feature = as.character(feature), 
    #         run = as.character(run), 
    #         label = as.character(label)
    #     )
    if ("censored" %in% names(processedData)) {
        df_allftr <- processedData %>%
            select(
                protein = PROTEIN,
                peptide = PEPTIDE,
                feature = FEATURE, 
                run = originalRUN, 
                label = LABEL,
                log2inty = ABUNDANCE, 
                is_censored = censored
            )
    } else {
        df_allftr <- processedData %>%
            select(
                protein = PROTEIN,
                peptide = PEPTIDE,
                feature = FEATURE, 
                run = originalRUN, 
                label = LABEL,
                log2inty = ABUNDANCE
            ) %>% 
            mutate(is_censored = FALSE)
    }
    df_allftr <- df_allftr %>% 
        mutate(
            is_obs = !(is.na(log2inty) | is_censored), 
            log2inty = ifelse(is_obs, log2inty, NA), 
            protein = as.character(protein), 
            peptide = as.character(peptide), 
            feature = as.character(feature), 
            run = as.character(run), 
            label = as.character(label)
        )
    
    # Is the dataset from a labeled experiment?
    is_labeled <- any(df_allftr$label == "H")
    all_label <- unique(df_allftr$label)
    
    # Process H/L channels separately
    df_rmftr <- vector("list", length(all_label))
    df_rmpk <- vector("list", length(all_label))
    for (l in seq_along(all_label)) {
        
        lab <- all_label[l]
        df_onelab <- df_allftr %>% filter(label == lab)
        message(paste0("Analyzing features in ", lab, " channel..."))
        
        # In heavy channel, only positive peaks are used for modeling
        if (lab == "H") {
            df_onelab <- df_onelab %>% filter(log2inty > 0)
        }
        
        # Singleton (unreplicated) or undetectable features
        ftr_unrep <- df_onelab %>% 
            group_by(protein, feature) %>% 
            filter(sum(is_obs) <= 1) %>% 
            ungroup() %>% 
            distinct(protein, feature)
        
        # Remove singleton or undetectable features & uncovered runs
        df_onelab <- df_onelab %>% 
            anti_join(ftr_unrep) %>% 
            group_by(protein, run) %>% 
            filter(any(is_obs)) %>% 
            ungroup()
        
        ### [TEST] ###
        # Proteins with sufficient features to be analyzed
        prot_nbftr <- df_onelab %>% 
            group_by(protein) %>% 
            summarise(nb_feature = n_distinct(feature)) %>% 
            filter(nb_feature > min_ftr)
        df_onelab <- df_onelab %>% 
            semi_join(prot_nbftr)
        ### [TEST] ###
        
        # Coverage of each feature
        obs_protein <- df_onelab %>% 
            group_by(protein) %>% 
            summarise(
                nb_run = n_distinct(run), 
                nb_feature = n_distinct(feature), 
                nb_obs = sum(is_obs)
            ) %>% 
            mutate(nb_full = nb_run * nb_feature, pi_obs = nb_obs / nb_full) %>% 
            mutate(min_obs = qbinom(0.01, nb_run, pi_obs))
        
        # Nested data frame
        nested_prot <- df_onelab %>% 
            select(protein, run, peptide, feature, log2inty, is_obs) %>% 
            group_by(protein) %>% 
            nest()
        
        nested_prot <- nested_prot %>% 
            left_join(obs_protein %>% select(protein, min_obs))
        
        # Coverage
        message("Identifying low-coverage features...")
        
        list_cover <- vector("list", length = nrow(nested_prot))
        for (i in seq_along(list_cover)) {
            augdata <- nested_prot$data[[i]] %>% 
                group_by(feature) %>% 
                mutate(is_lowcvr = sum(is_obs) < nested_prot$min_obs[i]) %>% 
                ungroup()
            nested_prot$data[[i]] <- augdata
            list_cover[[i]] <- augdata %>% distinct(feature, is_lowcvr)
        }
        nested_prot$cover_feature <- list_cover
        rm(list_cover)
        
        # List of low-coverage features
        ftr_lowcvr <- nested_prot %>% 
            unnest(cover_feature) %>% 
            filter(is_lowcvr) %>% 
            distinct(protein, feature)
        
        # Keep proteins with at least 3 good-coverage features
        ### [TEST] ###
        is_ftreff <- map_lgl(nested_prot$cover_feature, ~sum(!(.$is_lowcvr)) > min_ftr)
        # is_ftreff <- map_lgl(nested_prot$cover_feature, ~sum(!(.$is_lowcvr)) > 2)
        ### [TEST] ###
        nested_prot <- nested_prot[is_ftreff, ]
        
        # Robust linear model (NB: non-converged models are not used)
        message("Fitting robust linear models...")
        
        list_sdata <- map(nested_prot$data, ~filter(., !is_lowcvr))
        fit <- list_sdata %>% map(possibly(quietly(fit_prot_huber), NULL))
        is_failed <- map_lgl(fit, is.null)
        nested_prot <- nested_prot[!is_failed, ]
        fit <- fit[!is_failed] %>% transpose()
        is_warned <- !map_lgl(fit$warnings, is_empty)
        nested_prot$rlm_fit <- fit$result
        nested_prot <- nested_prot[!is_warned, ] 
        nested_prot <- nested_prot %>% 
            mutate(df_resid = map_dbl(rlm_fit, ~summary(.)$df[2])) %>% 
            mutate(s_resid = map_dbl(rlm_fit, ~summary(.)$sigma))
# nested_prot %>% 
#             slice(which(!(is_warned))) %>% 
#             mutate(df_resid = map_dbl(rlm_fit, ~summary(.)$df[2])) %>% 
#             mutate(s_resid = map_dbl(rlm_fit, ~summary(.)$sigma))
        rm(list = c("list_sdata", "fit", "is_failed", "is_warned"))
        
        # Extract variances and degrees of freedom
        var_protein <- nested_prot %>% 
            select(protein, s_resid, df_resid) %>% 
            mutate(var_resid = s_resid ^ 2) %>% 
            filter(!is.na(s_resid), df_resid > 0)
        
        # Shrinkage variance estimation with limma
        eb_fit <- limma::squeezeVar(var_protein$var_resid, var_protein$df_resid, robust = TRUE)
        var_protein <- var_protein %>% 
            ungroup() %>%
            mutate(var_resid_eb = eb_fit$var.post) %>% 
            mutate(s_resid_eb = sqrt(var_resid_eb))
        
        nested_prot <- nested_prot %>% 
            left_join(var_protein %>% select(protein, s_resid_eb))
        
        # Feature variance and outlier detection
        message("Identifying outliers and calculating feature variances...")
        
        list_var <- vector("list", nrow(nested_prot))
        for (i in seq_along(list_var)) {
            s_eb <- nested_prot$s_resid_eb[i]
            if (!is.na(s_eb)) {
                augdata <- nested_prot$rlm_fit[[i]] %>% augment()
                if (lab == "L") {
                    list_var[[i]] <- calc_fvar(augdata, s_eb, rm_olr = TRUE)
                    nested_prot$data[[i]] <- suppressMessages(
                        nested_prot$data[[i]] %>% 
                            left_join(flag_outlier(augdata, s_eb, keep_run = FALSE)) %>% 
                            mutate(is_olr = ifelse(is.na(is_olr), FALSE, is_olr))
                    )
                } else {
                    # No outlier detection in heavy channel; all measurements for calculation
                    list_var[[i]] <- calc_fvar(augdata, s_eb, rm_olr = FALSE)
                    nested_prot$data[[i]] <- nested_prot$data[[i]] %>% mutate(is_olr = FALSE)
                }
            } else {
                nested_prot$data[[i]] <- nested_prot$data[[i]] %>% mutate(is_olr = FALSE)
            }
        }
        nested_prot$var_feature <- list_var
        rm("list_var")
        
        # Noisy features
        allfvar <- nested_prot %>% 
            filter(!is.na(s_resid_eb)) %>% 
            select(protein, var_feature) %>% 
            unnest(var_feature)
        fvar_cut <- quantile(allfvar$svar_ref, 0.05, na.rm = TRUE)
        ftr_hivar <- allfvar %>% filter(svar_feature > fvar_cut) %>% distinct(protein, feature)
        
        # Combine unreplicated, low-coverage, and noisy features
        df_rmftr[[l]] <- bind_rows(ftr_hivar, ftr_lowcvr, ftr_unrep) %>% 
            select(protein, feature) %>% 
            mutate(label = lab)
        
        # Outliers
        df_rmpk[[l]] <- nested_prot %>% 
            unnest(data) %>% 
            filter(is_olr) %>% 
            select(protein, run, feature) %>% 
            mutate(label = lab)
        
        message(paste0("Outlier detection and feature selection in ", lab, " channel are completed"))
    }
    df_rmftr <- bind_rows(df_rmftr)
    df_rmpk <- bind_rows(df_rmpk)
    
    # Non-positive peaks in heavy channel are considered as outliers
    if (is_labeled) {
        df_rmpk <- bind_rows(
            df_rmpk, 
            df_allftr %>% 
                filter(label == "H", log2inty <= 0) %>% 
                select(protein, run, feature, label)
        )
    }
    
    # Annotate feature and peak qualities with 2 additional columns
    df_rmftr <- df_rmftr %>%
        rename(PROTEIN = protein, FEATURE = feature, LABEL = label)
    df_rmpk <- df_rmpk %>%
        rename(PROTEIN = protein, originalRUN = run, FEATURE = feature, LABEL = label)
    
    annotatedData <- bind_rows(
        suppressWarnings(processedData %>% semi_join(df_rmftr) %>% mutate(feature_quality = "Uninformative")),
        suppressWarnings(processedData %>% anti_join(df_rmftr) %>% mutate(feature_quality = "Informative"))
    )
    annotatedData <- bind_rows(
        suppressWarnings(annotatedData %>% semi_join(df_rmpk) %>% mutate(is_outlier = TRUE)),
        suppressWarnings(annotatedData %>% anti_join(df_rmpk) %>% mutate(is_outlier = FALSE))
    )
    
    return(annotatedData)
}
