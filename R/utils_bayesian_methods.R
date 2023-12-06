#' @importFrom protDP dpc
#' 
#' @export
MSstatsBayesSummarize = function(data, dpc_betas, 
                                 bayes_method="MCMC", 
                                 n_iterations=1000, 
                                 chains=4, cores=4, 
                                 elbo_samples=500,
                                 tol_rel_obj=.00001){
    
    keep = c("PROTEIN", "FEATURE", "RUN", "ABUNDANCE")
    subset_data = data[, ..keep]
    
    model_data = prepare_for_bayes(subset_data, dpc_betas)
    feature_data = model_data[[1]]
    stan_input = model_data[[2]]
    missing_runs = model_data[[3]]
    
    stan_file = system.file("stan", "bayes_model.stan", package="MSstats")
    
    if (bayes_method == "MCMC"){
        fit = stan(file = stan_file, data = stan_input,
                   chains = chains, iter = n_iterations, cores = cores,
                   seed=100, model_name="MSstats_model")
    } else if (bayes_method == "VB"){
        model = stan_model(stan_file)
        fit = vb(model, data = stan_input, iter = n_iterations, 
                 elbo_samples = elbo_samples, tol_rel_obj = tol_rel_obj)
    }
    model_results = recover_data(fit, feature_data, missing_runs)
    
    return(model_results)
}

calculate_dpc = function(data){
    
    wide_data = data.table::dcast(data, FEATURE~RUN, value.var = "ABUNDANCE")
    wide_data[, FEATURE := NULL]
    
    dpc_params = dpc(wide_data)
    
    return(dpc_params$beta)
}

MSstatsBayesSummarizationOutput = function(input, summarized){
    
    # Define feature and protein data
    feature_data = input
    protein_data = summarized
    
    # Add missing value info to feature level data
    feature_data = merge(feature_data, 
                         unique(protein_data[,c("PROTEIN_original", 
                                                "RUN_original",
                                                "FEATURE_original", 
                                                "imp_mean", "imp_sd")]),
                         by.x = c("PROTEIN", "RUN", "FEATURE"),
                         by.y = c("PROTEIN_original", "RUN_original", 
                                  "FEATURE_original"), 
                         all.x=TRUE, all.y=FALSE)
    feature_data[imp_mean == 0, imp_mean:=NA]
    feature_data[imp_sd == 0, imp_sd:=NA]
    
    # Calculate summary stats for protein level data
    summary_stats = feature_data
    summary_stats[, NonMissingStats := .getNonMissingFilterStats(.SD, NULL)]
    summary_stats[, NumMeasuredFeature := sum(NonMissingStats), 
          by = c("PROTEIN", "RUN")]
    summary_stats[, MissingPercentage := 1 - (NumMeasuredFeature / total_features)]
    summary_stats[, more50missing := MissingPercentage >= 0.5]
    summary_stats[, nonmissing_orig := LABEL == "L" & !censored]
    summary_stats[, NumImputedFeature := sum(LABEL == "L" & !nonmissing_orig),
          by = c("PROTEIN", "RUN")]
    summary_stats = summary_stats[, c("PROTEIN", "RUN", 
                                      "NumMeasuredFeature", "MissingPercentage",
                                      "more50missing", "NumImputedFeature")]
    
    
    # Format feature data
    feature_data[ ,c("GROUP", "SUBJECT", "nonmissing", "n_obs", 
                     "n_obs_run", "total_features", "prop_features") := NULL]
    
    setnames(feature_data, 
             c("GROUP_ORIGINAL", "SUBJECT_ORIGINAL", "imp_mean", "imp_sd"), 
             c("GROUP", "SUBJECT", "predicted", "predicted_sd"))
    
    feature_data = feature_data[, c("PROTEIN", "PEPTIDE", "TRANSITION", 
                                    "FEATURE", "LABEL", "GROUP", "RUN", 
                                    "SUBJECT", "FRACTION", "originalRUN", 
                                    "censored", "INTENSITY", "ABUNDANCE",
                                    "newABUNDANCE", "predicted", "predicted_sd"
                                    )]
    
    
    # Add run level summarization
    summary_cols = c("PROTEIN_original", "RUN_original", "run_mean", "run_sd")
    protein_data = merge(unique(protein_data[, ..summary_cols]), 
                         unique(input[, c("PROTEIN", "RUN", "originalRUN", 
                                          "GROUP_ORIGINAL", "SUBJECT_ORIGINAL")]
                                ),
                         by.x = c("PROTEIN_original", "RUN_original"), 
                         by.y = c("PROTEIN", "RUN"), 
                         all.x=TRUE, all.y=FALSE)
    
    protein_data = merge(protein_data, unique(summary_stats),
                         by.x = c("PROTEIN_original", "RUN_original"),
                         by.y = c("PROTEIN", "RUN"), 
                         all.x=TRUE, all.y=FALSE)
    
    
    protein_data = merge(protein_data, 
                         feature_data[, uniqueN(.SD),
                                      by = c("PROTEIN", "GROUP"),
                                      .SDcols = c("FEATURE", "originalRUN")],
                         by.x = c("PROTEIN_original", "GROUP_ORIGINAL"),
                         by.y=c("PROTEIN", "GROUP"),
                         all.x=TRUE, all.y=FALSE)
    
    
    setnames(protein_data, 
             c("PROTEIN_original", "GROUP_ORIGINAL", "RUN_original", 
               "run_mean", "run_sd", "SUBJECT_ORIGINAL", "V1"), 
             c("Protein", "GROUP", "RUN", "LogIntensities", 
               "LogIntensities_sd", "SUBJECT", "TotalGroupMeasurements"))
    
    protein_data = protein_data[,c("RUN", "Protein", "LogIntensities", 
                                   "LogIntensities_sd", "originalRUN", "GROUP",
                                   "SUBJECT", "TotalGroupMeasurements", 
                                   "NumMeasuredFeature", "MissingPercentage", 
                                   "more50missing", "NumImputedFeature")]
    
    
    return(list(FeatureLevelData = as.data.frame(feature_data), 
               ProteinLevelData = as.data.frame(protein_data), 
               SummaryMethod = "Bayesian")
           )
}

prepare_for_bayes = function(data, dpc_betas){
    
    data_list = remove_missing_runs(data)
    data = data_list[[1]]
    missing_data = data_list[[2]]
    
    format_data = format_data(data)
    priors = get_priors(data)
    results = flatten_input(format_data)
    
    input = results[[1]]
    lookup_table = results[[2]]
    
    # grab summarized priors
    run_priors = merge(priors,
                       lookup_table[,c("PROTEIN_original", 
                                       "RUN_original",
                                       "Protein_run_idx")],
                       all.x=TRUE, all.y=FALSE, 
                       by=c("PROTEIN_original", "RUN_original"))
    run_priors = arrange(unique(run_priors), Protein_run_idx)
    
    stan_input = arrange_stan_data(input, run_priors, dpc_betas)
    
    return(list(lookup_table, stan_input, missing_data))
}

remove_missing_runs = function(data){
    
    data[, run_total:=all(is.na(ABUNDANCE)), by=.(PROTEIN, RUN)]
    
    remove = data[run_total != 0,][, run_total := NULL]
    keep = data[run_total == 0,][, run_total := NULL]
    
    return(list(keep, remove))
    
}

arrange_stan_data = function(input, run_priors, dpc_betas){
    num_proteins = length(unique(input$PROTEIN))
    num_runs = length(unique(input$Protein_run_idx))
    num_feat = length(unique(input$Protein_feature_idx))
    
    obs_mat = matrix(as.numeric(unlist(input[which(!is.na(input[,3])),])),
                     nrow=nrow(input[which(!is.na(input[,3])),]))
    missing_mat = matrix(as.numeric(unlist(input[which(is.na(input[,3])),])),
                         nrow=nrow(input[which(is.na(input[,3])),]))
    
    stan_input = list(N_obs=nrow(obs_mat), 
                      N_missing=nrow(missing_mat),
                      R=num_runs,
                      Feat=num_feat,
                      P=num_proteins,
                      obs=obs_mat[,3],
                      run_id=obs_mat[,1],
                      feature_id=obs_mat[,2],
                      protein_id=obs_mat[,4],
                      run_id_missing=missing_mat[,1],
                      feature_id_missing=missing_mat[,2],
                      protein_id_missing=missing_mat[,4],
                      zeros=rep(0, nrow(obs_mat)),
                      ones=rep(0, nrow(missing_mat)),
                      run_mu_prior=run_priors[,3][[1]],
                      sigma_run_prior=run_priors[,4][[1]],
                      feat_mu_prior=rep(0,num_feat),
                      sigma_feat_prior=rep(1, num_feat),
                      sigma_prior=rep(2, num_runs))#,
                      # beta0=dpc_betas[[1]],
                      # beta1=dpc_betas[[2]])
}

format_data = function(data){
    
    # Missing indicator
    data[, "Missing"] = ifelse(is.na(data[, "ABUNDANCE"]), 1., 0.)
    
    for (col in c("PROTEIN", "RUN", "FEATURE")){
        id_map = as.numeric(droplevels(unlist(data[, ..col])))
        data[, paste0(col, "_original")] = data[,..col]
        data[,col] = id_map
    }
    
    # order
    data = data[with(data, order(PROTEIN, FEATURE, RUN)), ]
    
    return(data)
}

flatten_input = function(data) {
    
    # Convert the data to a data frame
    # data = data.frame(data, stringsAsFactors = FALSE)
    # colnames(data) = c("Protein", "Run", "Feature", "Intensity", "Missing")
    
    # Add a list index
    # data$list_index = 1:nrow(data)
    
    # Create a unique protein_run identifier
    data$Protein_run = paste(data$PROTEIN, data$RUN, sep = "_")
    
    # Create a unique protein_feature identifier
    data$Protein_feature = paste(data$PROTEIN, data$FEATURE, sep = "_")
    
    # Create protein_run_idx and protein_feature_idx
    protein_run_idx = unique(data$Protein_run)
    protein_feature_idx = unique(data$Protein_feature)
    
    data$Protein_run_idx = match(data$Protein_run, protein_run_idx)
    data$Protein_feature_idx = match(data$Protein_feature, protein_feature_idx)
    
    # Return the flattened data
    flatten_data = data[, c("Protein_run_idx", "Protein_feature_idx", 
                            "ABUNDANCE", "PROTEIN")]
    
    return(list(flatten_data, data))
}

get_priors = function(data) {
    
    # Initialize collection structures
    run_priors = list()
    run_priors_std = list()
    feature_priors = list()
    
    proteins = unique(data[, 1])[[1]]
    
    for (i in seq_along(proteins)) {
        prot = as.character(proteins[i])
        temp_data = data[which(data[, 1] == prot)]
        
        runs = unique(temp_data[, 3])[[1]]
        # features = unique(temp_data[, 8])
        
        # Overall mean
        overall_mean = mean(temp_data[, 4][[1]], na.rm=TRUE)
        
        # Calculate run priors
        run_effect = c()
        run_std = list()
        
        for (r in seq_along(runs)) {
            run_effect = c(run_effect, 
                   mean(temp_data[
                   which(temp_data[, 3] == as.numeric(runs[r])), 4][[1]], 
                   na.rm=TRUE))
            run_std = c(run_std,
                        sd(temp_data[
                            which(temp_data[, 3] == as.numeric(runs[r])), 4][[1]], 
                            na.rm=TRUE))
        }
        
        run_effect = as.numeric(run_effect)
        run_std = as.numeric(run_std)
        run_std[is.na(run_std)] = 2.
        temp_df = data.frame("PROTEIN_original" = rep(prot, length(runs)),
                             "RUN_original" = runs,
                             "Run_prior" = run_effect,
                             "Run_std" = run_std
                             )
        # run_effect[is.na(run_effect)] = 0
        run_priors[[i]] = temp_df[which(!is.na(temp_df[,3])),]
        
        # run_std = as.numeric(run_std)
        # run_std[is.na(run_std)] = 0.1
        # run_std[run_std == 0] = 0.75
        # run_priors_std = c(run_priors_std, run_std)
        
        # # Calculate feature priors
        # feature_effect = list()
        # feature_std = list()
        # 
        # for (f in features) {
        #   feature_effect = c(feature_effect,
        #                       mean(temp_data[, 3][(temp_data[, 2] == f) & (temp_data[, 3] != 0)] - overall_mean))
        #   feature_std = c(feature_std, sd(temp_data[, 3][(temp_data[, 2] == f) & (temp_data[, 3] != 0)]))
        # }
        # 
        # feature_effect = as.numeric(feature_effect)
        # feature_effect[is.na(feature_effect)] = 0
        # feature_priors = c(feature_priors, feature_effect)
    }
    
    # priors = list(run_priors = unlist(run_priors),
    #                run_priors_std = unlist(run_priors_std))
    # ,feature_priors = feature_priors)
    
    return(data.table::rbindlist(run_priors))
}

recover_data = function(fit, feature_data, missing_runs){
    
    parameters = c("obs_mis", "run_mu", "sigma", "beta0", "beta1")
    
    # Missing imputation
    results = rstan::summary(fit, pars="obs_mis")
    results = results$summary[,c("mean", "sd")]
    imp_index = which(is.na(feature_data$ABUNDANCE))
    
    feature_data[, "imp_mean"] = 0
    feature_data[, "imp_sd"] = 0
    feature_data[imp_index, "imp_mean" := results[,1]]
    feature_data[imp_index, "imp_sd" := results[,2]]
    
    # Feature estimation
    results = rstan::summary(fit, pars="feature_mu")
    results = as.data.frame(results$summary[,c("mean", "sd")])
    results = rename(results, c("feature_mean" = "mean", 
                                "feature_sd" = "sd"))
    results$Protein_feature_idx = seq_len(nrow(results))
    feature_data = merge(feature_data, results, all.x=TRUE, 
                         by="Protein_feature_idx")
    
    # Run estimation
    results = rstan::summary(fit, pars="run_mu")
    results = as.data.frame(results$summary[,c("mean", "sd")])
    results = rename(results, c("run_mean" = "mean", 
                                "run_sd" = "sd"))
    results$Protein_run_idx = seq_len(nrow(results))
    feature_data = merge(feature_data, results, all.x=TRUE, 
                         by="Protein_run_idx")
    
    
    # beta0 = rstan::summary(fit, pars="beta0")
    # beta1 = rstan::summary(fit, pars="beta1")
    run_mu = rstan::summary(fit, pars="run_mu")
    feature_mu = rstan::summary(fit, pars="feature_mu")
    obs_mis = rstan::summary(fit, pars="obs_mis")
    sigma = rstan::summary(fit, pars="sigma")
    
    return(list(result_df = feature_data,
                bayes_results=list(
                    # beta = rbind(beta0$summary, beta1$summary),
                    run = run_mu$summary,
                    feature = feature_mu$summary,
                    missing = obs_mis$summary,
                    sigma = sigma$summary)
                )
           )
}