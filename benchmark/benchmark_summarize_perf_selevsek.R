library(MSstats)

# Performance expectations for MSstatsSummarizeWithMultipleCores on this dataset
# (not enforced automatically here -- read the job output and compare by hand):
#   - 4 cores should reduce wall time by at least 25% vs. 1 core
#   - peak RSS per worker should stay under ~1GB (see the memory report below)

input <- data.table::fread("/projects/VitekLab/Data/MS/selevsek/before_summarization.csv")

cat("=== MSstatsSummarizeWithMultipleCores performance check ===\n")
cat("Expectation: >=25% reduction in wall time on 4 cores vs. 1 core.\n")
cat("Expectation: peak RSS per worker should stay under ~1GB (see memory report below).\n\n")

time_1core <- system.time(
  result_1core <- MSstatsSummarizeWithSingleCore(
    input, "TMP", TRUE, "NA", FALSE, TRUE
  )
)

time_4core <- system.time(
  result_4core <- MSstatsSummarizeWithMultipleCores(
    input, "TMP", TRUE, "NA", FALSE, TRUE, 4, track_memory = TRUE,
    max_proteins_per_worker = 200
  )
)

elapsed_1core <- time_1core[["elapsed"]]
elapsed_4core <- time_4core[["elapsed"]]
speedup_pct <- (elapsed_1core - elapsed_4core) / elapsed_1core * 100

cat(sprintf("1 core wall time:  %.2f s\n", elapsed_1core))
cat(sprintf("4 core wall time:  %.2f s\n", elapsed_4core))
cat(sprintf("Observed speedup:  %.1f%% (expectation: >= 25%%)\n", speedup_pct))
cat("Compare the 'Worker RSS min / mean / max' line above against the ~1GB expectation.\n")

cat("\n=== Verifying 1-core and 4-core results are identical ===\n")

compare_df <- function(df1, df2) {
  df1 <- as.data.frame(df1)
  df2 <- as.data.frame(df2)

  if (!setequal(names(df1), names(df2))) return(FALSE)

  df2 <- df2[, names(df1), drop = FALSE]

  # Convert factors to character so ordering/comparison is based on
  # content, not on (possibly differing) factor level order
  df1[] <- lapply(df1, function(x) if (is.factor(x)) as.character(x) else x)
  df2[] <- lapply(df2, function(x) if (is.factor(x)) as.character(x) else x)

  ord1 <- do.call(order, as.list(df1))
  ord2 <- do.call(order, as.list(df2))

  df1 <- df1[ord1, , drop = FALSE]
  df2 <- df2[ord2, , drop = FALSE]

  rownames(df1) <- NULL
  rownames(df2) <- NULL

  isTRUE(all.equal(df1, df2, tolerance = 1e-14))
}

n_proteins <- length(result_1core)

protein_level_matches <- logical(n_proteins)
for (i in seq_len(n_proteins)) {
  protein_level_matches[i] <- compare_df(result_1core[[i]][[1]], result_4core[[i]][[1]])
}

feature_level_data_matches <- logical(n_proteins)
for (i in seq_len(n_proteins)) {
  feature_level_data_matches[i] <- compare_df(result_1core[[i]][[2]], result_4core[[i]][[2]])
}

n_protein_level_match <- sum(protein_level_matches)
n_feature_level_data_match <- sum(feature_level_data_matches)

cat(sprintf("Protein-level results matching: %d / %d\n", n_protein_level_match, n_proteins))
cat(sprintf("Feature_level data matching:         %d / %d\n", n_feature_level_data_match, n_proteins))

stopifnot(
  "Protein-level results differ between 1 core and 4 cores" =
    n_protein_level_match == n_proteins,
  "Feature_level data differs between 1 core and 4 cores" =
    n_feature_level_data_match == n_proteins
)

cat("1-core and 4-core results are identical.\n")
