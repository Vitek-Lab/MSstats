#!/usr/bin/env Rscript
#
# profile_dataprocess_peak.R
#
# Stage-by-stage peak-memory instrumentation for MSstats::dataProcess.
#
# Replays dataProcess's internals manually with a checkpoint between each
# stage and prints the high-water mark (gc()'s "max used" column) at
# every boundary. Useful for locating which pipeline stage drives overall
# peak memory -- complementary to profmem::profmem(), which pinpoints
# the biggest single allocation but can miss cases where peak is caused
# by many small coexisting objects.
#
# Usage
# -----
#
#   # From the package root (where DESCRIPTION lives):
#   Rscript benchmark/profile_dataprocess_peak.R
#
#   # With a larger fixture -- the argument is the DDARawData replication
#   # factor (default 100):
#   Rscript benchmark/profile_dataprocess_peak.R 300
#
# Requirements
# ------------
#
#   - pkgload must be installed.
#   - Working directory must be the package root so pkgload::load_all(".")
#     can find DESCRIPTION and load the branch under test.
#
# How to read the output
# ----------------------
#
# Two columns per row:
#   used : Vcells currently in use at the checkpoint (can shrink when
#          intermediates are freed and gc() runs).
#   max  : running maximum of `used` since gc(reset = TRUE) at "start".
#          Never decreases.
#
# The *delta in the max column between consecutive rows* is the peak
# contribution of that stage. Zero delta = the stage allocated but
# stayed below the prior high-water mark. Positive delta = the stage
# bumped the mark.
#
# Example output (n_replicates=100, ~10.8 MB input):
#
#   start                                     used= 47.4  max= 47.4
#   after PrepareForDataProcess + rm(raw)     used= 61.3  max=155.4   <- +108 MB
#   after Normalize                           used= 61.3  max=155.4
#   after MergeFractions                      used= 61.3  max=155.4
#   after HandleMissing                       used= 62.1  max=155.4
#   after SelectFeatures                      used= 62.2  max=155.4
#   after PrepareForSummarization             used= 73.2  max=155.4
#   after Summarize                           used=102.1  max=155.6
#   after gc()                                used=102.1  max=155.6
#   after SummarizationOutput                 used= 92.4  max=178.9   <- +23 MB
#
# Interpretation: the ~131 MB peak-above-baseline lives in two transient
# spikes -- ~108 MB in MSstatsPrepareForDataProcess at pipeline entry,
# and ~23 MB more in MSstatsSummarizationOutput at exit. Stages in
# between never breach the mark Prepare set.
#
# Caveats
# -------
#
#   - Uses MSstats::: (triple colon) to reach unexported helpers. Only
#     OK because this script is for diagnosis, not for shipping.
#   - The stage sequence mirrors R/dataProcess.R and must be kept in
#     sync if that pipeline is reordered.
#   - Run-to-run noise of about +/- 10 MB is normal (GC timing, OS page
#     caches). The relative deltas between stages are stable.
#
# -----------------------------------------------------------------------

# --- CLI arg: number of DDARawData replicates (fixture size) ------------

args <- commandArgs(trailingOnly = TRUE)
n_replicates <- if (length(args) >= 1) as.integer(args[1]) else 100L
if (is.na(n_replicates) || n_replicates < 1) {
    stop("Expected a positive integer n_replicates as the first argument")
}

# --- Load the package from source ---------------------------------------

if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required: install.packages('pkgload')")
}
pkgload::load_all(".", quiet = TRUE)
MSstatsConvert::MSstatsLogsSettings(FALSE)

# --- Build the fixture --------------------------------------------------
#
# Same recipe the memory tests use: replicate the built-in DDARawData and
# suffix ProteinName so each replicate is treated as a distinct protein.

set.seed(42)
base_data <- data.table::as.data.table(DDARawData)
replicated_data <- data.table::rbindlist(lapply(seq_len(n_replicates), function(i) {
    d <- data.table::copy(base_data)
    d$ProteinName <- paste0(d$ProteinName, "_rep", i)
    d
}))
replicated_data <- as.data.frame(replicated_data)
input_mb <- as.numeric(object.size(replicated_data)) / 1e6

cat(sprintf("Fixture: %d rows, %.1f MB  (n_replicates=%d)\n\n",
            nrow(replicated_data), input_mb, n_replicates))

# --- Checkpoint helper --------------------------------------------------
#
# gc()'s matrix layout varies across R versions (some add a "limit (Mb)"
# column in 4.x). Find the "(Mb)" column that immediately follows each
# labelled count column by name rather than by hard-coded position.

checkpoint <- function(label, reset = FALSE) {
    if (reset) { gc(); gc(reset = TRUE) }
    g <- gc()
    cols <- colnames(g)
    used_mb <- g["Vcells", which(cols == "used")[1] + 1]
    max_mb  <- g["Vcells", which(cols == "max used")[1] + 1]
    cat(sprintf("%-40s  used=%6.1f  max=%6.1f\n", label, used_mb, max_mb))
}

# --- Replay dataProcess stages manually ---------------------------------
#
# Mirrors the sequence in R/dataProcess.R::dataProcess. Keep these in
# sync if dataProcess reorders its pipeline.

checkpoint("start", reset = TRUE)

raw <- replicated_data
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
