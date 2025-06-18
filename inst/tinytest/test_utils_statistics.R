# Test .logMissingness
logged_msgs <- list()
options(MSstatsLog = function(level, msg) logged_msgs <<- c(logged_msgs, msg))
options(MSstatsMsg = function(level, msg) {})  # no-op

input <- data.table::data.table(
  LABEL = rep(c("L", "H"), each = 12),
  GROUP = rep(c("Control", "Treated"), each = 6, times = 2),
  FEATURE = rep(c("FEATURE1", "FEATURE2", "FEATURE3"), times = 8),
  RUN = rep(paste0("Run", 1:4), each = 3, times = 2),
  ABUNDANCE = c(
    # Control group
    10.5, 11.2,  NA,     # Run1
    10.8,  NA,   NA,     # Run2 (2 missing)
    11.0, 11.1, 11.2,    # Run3
    NA,   NA,   NA,      # Run4 (fully missing)
    # Treatment group
    9.5, 9.8,   10.0,    # Run1
    NA,   NA,   NA,      # Run2 (fully missing)
    9.2, 9.3,    NA,     # Run3
    NA,   NA,   NA       # Run4 (fully missing)
  )
)

MSstats:::.logMissingness(input)

expect_true(
  any(grepl("Some features are completely missing.*PEPTIDE2|PEPTIDE3", logged_msgs)),
  info = "Should log that at least one feature is completely missing in a condition"
)

expect_true(
  any(grepl("more than 75% missing values.*Run4", logged_msgs)) &&
    any(grepl("Run2", logged_msgs)),
  info = "Should log that Run2 and Run4 have >75% missing"
)
