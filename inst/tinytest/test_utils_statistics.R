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