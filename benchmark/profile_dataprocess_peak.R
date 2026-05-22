#!/usr/bin/env Rscript
#
# profile_dataprocess_peak.R
#
# Peak-memory + wall-clock instrumentation for MSstats::dataProcess.
#
# Replays dataProcess's internals manually, printing the gc() high-water
# mark and elapsed time at every stage boundary. Useful for locating
# which pipeline stage drives overall peak memory and runtime.
#
# Usage
# -----
#
#   # From the package root:
#   Rscript benchmark/profile_dataprocess_peak.R
#
#   # With an alternate MSstats-format CSV (10 cols: ProteinName,
#   # PeptideSequence, PrecursorCharge, FragmentIon, ProductCharge,
#   # IsotopeLabelType, Run, BioReplicate, Condition, Intensity):
#   Rscript benchmark/profile_dataprocess_peak.R path/to/input.csv
#
# How to read the output
# ----------------------
#
#   used : Vcells currently in use at the checkpoint.
#   max  : running max of `used` since gc(reset = TRUE) at "start".
#   dt   : seconds elapsed since the previous checkpoint.
#
# The *delta in the max column between consecutive rows* is the peak
# memory contribution of that stage. The dt column is per-stage wall time.

# --- CLI arg: path to MSstats-format CSV, or integer for DDARawData ----
#
# - Rscript profile_dataprocess_peak.R                 -> default CSV
# - Rscript profile_dataprocess_peak.R path/to/x.csv   -> use that CSV
# - Rscript profile_dataprocess_peak.R 100             -> replicate
#                                                         DDARawData 100x

default_input <- "data/reduce_output_spectronaut_input_qc_test.csv"
args <- commandArgs(trailingOnly = TRUE)
arg <- if (length(args) >= 1) args[1] else default_input
n_replicates <- suppressWarnings(as.integer(arg))
use_ddaraw <- !is.na(n_replicates) && n_replicates >= 1
if (!use_ddaraw && !file.exists(arg)) {
    stop(sprintf("Input file not found: %s", arg))
}

# --- Load the package from source ---------------------------------------

if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required: install.packages('pkgload')")
}
pkgload::load_all(".", quiet = TRUE)
MSstatsConvert::MSstatsLogsSettings(FALSE)

# --- Load the fixture ---------------------------------------------------

required_cols <- c("ProteinName", "PeptideSequence", "PrecursorCharge",
                   "FragmentIon", "ProductCharge", "IsotopeLabelType",
                   "Run", "BioReplicate", "Condition", "Intensity")
if (use_ddaraw) {
    set.seed(42)
    base <- data.table::as.data.table(DDARawData)
    raw_input <- data.table::rbindlist(lapply(seq_len(n_replicates), function(i) {
        d <- data.table::copy(base)
        d$ProteinName <- paste0(d$ProteinName, "_rep", i)
        d
    }))
    raw_input <- as.data.frame(raw_input)
    fixture_label <- sprintf("DDARawData x %d", n_replicates)
} else {
    raw_input <- data.table::fread(arg, data.table = FALSE)
    fixture_label <- arg
    # If this looks like a raw Spectronaut export (has R.FileName /
    # PG.ProteinAccessions, etc.), convert to MSstats 10-col format.
    if ("R.FileName" %in% colnames(raw_input) &&
        !all(required_cols %in% colnames(raw_input))) {
        cat("Detected raw Spectronaut export — converting via",
            "SpectronauttoMSstatsFormat()...\n")
        raw_input <- MSstatsConvert::SpectronauttoMSstatsFormat(
            raw_input, use_log_file = FALSE)
        fixture_label <- paste0(arg, " (Spectronaut-converted)")
    }
}
missing_cols <- setdiff(required_cols, colnames(raw_input))
if (length(missing_cols) > 0) {
    stop(sprintf("Input is missing MSstats columns: %s",
                 paste(missing_cols, collapse = ", ")))
}
input_mb <- as.numeric(object.size(raw_input)) / 1e6
cat(sprintf("Fixture: %s\n         %d rows, %.1f MB in memory\n\n",
            fixture_label, nrow(raw_input), input_mb))

# --- Checkpoint helper --------------------------------------------------

bench_rows <- list()
.last_time <- NULL
checkpoint <- function(label, reset = FALSE) {
    if (reset) { gc(); gc(reset = TRUE); .last_time <<- proc.time() }
    now <- proc.time()
    dt  <- if (is.null(.last_time)) 0 else (now - .last_time)[["elapsed"]]
    .last_time <<- now
    g <- gc()
    cols <- colnames(g)
    used_mb <- g["Vcells", which(cols == "used")[1] + 1]
    max_mb  <- g["Vcells", which(cols == "max used")[1] + 1]
    cat(sprintf("%-42s  used=%7.1f  max=%7.1f  dt=%6.2fs\n",
                label, used_mb, max_mb, dt),
        file = stderr())
    bench_rows[[length(bench_rows) + 1L]] <<- data.frame(
        stage   = label,
        used_mb = used_mb,
        max_mb  = max_mb,
        dt_s    = dt,
        stringsAsFactors = FALSE
    )
    invisible(dt)
}

# --- Replay dataProcess stages manually ---------------------------------
#
# Mirrors the sequence in R/dataProcess.R::dataProcess. Keep these in
# sync if dataProcess reorders its pipeline.

t_start <- proc.time()
checkpoint("start", reset = TRUE)

raw <- raw_input
peptides_dict <- MSstats:::makePeptidesDictionary(
    data.table::as.data.table(unclass(raw)), "equalizeMedians")
input <- MSstats:::MSstatsPrepareForDataProcess(raw, 2, NULL)
rm(raw); checkpoint("after PrepareForDataProcess + rm(raw)")

input <- MSstats:::MSstatsNormalize(input, "equalizeMedians", peptides_dict, NULL)
rm(peptides_dict); checkpoint("after Normalize")

input <- MSstats:::MSstatsMergeFractions(input)
checkpoint("after MergeFractions")

input <- MSstats:::MSstatsHandleMissing(input, "TMP", TRUE, "NA", 0.999)
checkpoint("after HandleMissing")

input <- MSstats:::MSstatsSelectFeatures(input, "all", 3, 2)
checkpoint("after SelectFeatures")

processed <- MSstats:::getProcessed(input)
input <- MSstats:::MSstatsPrepareForSummarization(input, "TMP", TRUE, "NA", TRUE)
checkpoint("after PrepareForSummarization")

summarized <- MSstats:::MSstatsSummarizeWithMultipleCores(
    input, "TMP", TRUE, "NA", FALSE, TRUE, 1, 90)
checkpoint("after Summarize")

gc(verbose = FALSE)
checkpoint("after gc()")

output <- MSstats:::MSstatsSummarizationOutput(
    input, summarized, processed, "TMP", TRUE, "NA")
checkpoint("after SummarizationOutput")

total_elapsed <- (proc.time() - t_start)[["elapsed"]]

# --- Summary ------------------------------------------------------------

bench_df      <- do.call(rbind, bench_rows)
peak_mb       <- max(bench_df$max_mb)
baseline_mb   <- bench_df$max_mb[bench_df$stage == "start"]
peak_delta_mb <- peak_mb - baseline_mb

cat("\n========================================\n")
cat(sprintf("Total elapsed:        %.2f s\n",  total_elapsed))
cat(sprintf("Peak memory (max):    %.1f MB\n", peak_mb))
cat(sprintf("Baseline memory:      %.1f MB\n", baseline_mb))
cat(sprintf("Peak above baseline:  %.1f MB\n", peak_delta_mb))
cat("========================================\n")
