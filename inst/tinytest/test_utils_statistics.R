# Test .logSummaryStatistics
log_env <- new.env()
log_env$messages <- character()

options(MSstatsLog = function(level, msg) {
  log_env$messages <- c(log_env$messages, msg)
})
options(MSstatsMsg = function(level, msg) {})  # no-op for this test

## Simulate a small, realistic MS data input
input <- data.table::data.table(
  PROTEIN = c("P1", "P1", "P1", "P2", "P2", "P3"),
  PEPTIDE = c("pep1", "pep1", "pep2", "pep3", "pep3", "pep4"),
  FEATURE = c("f1", "f2", "f3", "f4", "f5", "f6"),
  RUN = c("Run1", "Run2", "Run3", "Run1", "Run2", "Run1"),
  SUBJECT_ORIGINAL = c("S1", "S1", "S2", "S3", "S3", "S4"),
  FRACTION = c(1, 1, 1, 1, 1, 1),
  GROUP_ORIGINAL = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")
)

## Run the function
MSstats:::.logSummaryStatistics(input)

## Extract logged messages
msgs <- log_env$messages

## Check for summary stats
expect_true(
  any(grepl("# proteins:\\s*3", msgs)),
  info = "Should log the correct number of proteins (3)"
)
expect_true(
  any(grepl("# peptides per protein:\\s*1-2", msgs)),
  info = "Should log the peptide per protein range (1–2)"
)
expect_true(
  any(grepl("# features per peptide:\\s*1-2", msgs)),
  info = "Should log the feature per peptide range (1–2)"
)

## Check for message about proteins with only one feature (P3)
expect_true(
  any(grepl("Some proteins have only one feature.*P3", msgs)),
  info = "Should log that P3 has only one feature"
)

## Check sample-level summary
expect_true(
  any(grepl("# runs", msgs)) &&
    any(grepl("# bioreplicates", msgs)) &&
    any(grepl("# tech\\. replicates", msgs)),
  info = "Should include runs, bio replicates, and tech replicates summary"
)

# Test .logMissingness
log_env <- new.env()
log_env$messages <- character()

## Override logging options with closures that update the environment
options(MSstatsLog = function(level, msg) {
    log_env$messages <- c(log_env$messages, msg)
})
options(MSstatsMsg = function(level, msg) {})  # optional no-op

## Realistic test dataset
input <- data.table::data.table(
    LABEL = rep(c("L", "H"), each = 12),
    GROUP = rep(c("Control", "Treatment"), each = 6, times = 2),
    FEATURE = rep(c("PEPTIDE1", "PEPTIDE2", "PEPTIDE3"), times = 8),
    RUN = rep(paste0("Run", 1:4), each = 3, times = 2),
    ABUNDANCE = c(
        10.5, 11.2,  NA,     # Run1
        10.8,  NA,   NA,     # Run2
        11.0, 11.1, 11.2,    # Run3
        NA,   NA,   NA,      # Run4
        9.5, 9.8,   10.0,    # Run1
        NA,   NA,   NA,      # Run2
        9.2, 9.3,    NA,     # Run3
        NA,   NA,   NA       # Run4
    )
)

## Run the function
MSstats:::.logMissingness(input)

## Extract messages
msgs <- log_env$messages

## Feature-level logging check
expect_true(
    any(grepl("Some features are completely missing.*PEPTIDE", msgs)),
    info = "Should log about completely missing features"
)

## Run-level logging check
expect_true(
    any(grepl("more than 75% missing values.*Run2", msgs)) &&
        any(grepl("Run4", msgs)),
    info = "Should log that Run2 and Run4 have >75% missing"
)