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
