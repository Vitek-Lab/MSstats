
.normalize_anomalies = function(data){
    # Step 1: Find the global minimum anomaly score
    global_min = min(data$ANOMALYSCORES )
    
    # Step 2: Compute the local minimum for each PSM and adjust scores
    data[, local_min := min(ANOMALYSCORES ), by = PEPTIDE]
    data[, normalized_score := ANOMALYSCORES - local_min + global_min]
    
    return(data)
}

.calculate_weights = function(data){
    
    data = .normalize_anomalies(data)

    w_stats = data[, .(
        shapiro_W = tryCatch(
            shapiro.test(unlist(normalized_score))$statistic,
            error = function(e) NA_real_
        )
    ), by = PEPTIDE]
    
    # Join the W stats back into the original data.table
    data = merge(data, w_stats, by = "PEPTIDE", all.x = TRUE)
    
    anomaly_scores = data$normalized_score
    imputation_var = data$imputation_var
    norm_score = data$shapiro_W
    
    if (sum(is.na(imputation_var)) == length(imputation_var)){
        weights = anomaly_scores
    } else {
        resPPCA = pca(matrix(
            c(anomaly_scores, imputation_var, norm_score), ncol=3), 
            method="ppca", center=TRUE, nPcs=1,
            seed=1, maxIterations = 10000)
        weights = (abs(resPPCA@loadings[[1]]) * anomaly_scores
        ) + 
            (abs(resPPCA@loadings[[3]]) * norm_score
            ) + ifelse(
                is.na(imputation_var), 0, 
                abs(resPPCA@loadings[[2]]) * imputation_var)
    }
    
    data$raw_weights = weights
    data$weights = 1 / weights
    
    return(data)
}