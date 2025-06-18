# Test .logMissingness
logged_msgs <- list()
options(MSstatsLog = function(level, msg) logged_msgs <<- c(logged_msgs, msg))
options(MSstatsMsg = function(level, msg) {})  # no-op

input <- data.table::data.table(
  LABEL = rep(c("A", "B"), each = 4),
  GROUP = rep(c("G1", "G2"), each = 2, times = 2),
  FEATURE = c(rep("F1", 4), rep("F2", 4)),
  ABUNDANCE = c(NA, NA, 1.2, 1.4, NA, NA, NA, NA),  # F2 is fully NA in all rows
  RUN = paste0("R", 1:8)
)

.logMissingness(input)

expect_true(
  any(grepl("Some features are completely missing.*F2", logged_msgs)),
  info = "Should log that F2 is completely missing in at least one condition"
)

expect_true(
  any(grepl("more than 75% missing values.*R5", logged_msgs)) &&
    any(grepl("R6", logged_msgs)) &&
    any(grepl("R7", logged_msgs)) &&
    any(grepl("R8", logged_msgs)),
  info = "Should log that runs R5 to R8 have >75% missing"
)
