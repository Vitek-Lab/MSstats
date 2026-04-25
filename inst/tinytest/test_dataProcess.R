# Test dataProcess with default parameters ------------------------------------
QuantDataDefault = dataProcess(SRMRawData, use_log_file = FALSE)
QuantDataDefaultLinear = dataProcess(DDARawData, use_log_file = FALSE,
                                     summaryMethod = "linear")

# Test dataProcess with numberOfCores parameter ----------------------
QuantDataParallel = dataProcess(SRMRawData, use_log_file = FALSE,
                                               numberOfCores = 2)
QuantDataParallelLinear = dataProcess(DDARawData, use_log_file = FALSE,
                                     summaryMethod = "linear", numberOfCores = 2)

expect_equal(nrow(QuantDataDefault$FeatureLevelData),
             nrow(QuantDataParallel$FeatureLevelData))

expect_equal(nrow(QuantDataDefaultLinear$FeatureLevelData),
             nrow(QuantDataParallelLinear$FeatureLevelData))


# Test multicore linear parity --------------------------------------------
# Linear summarization is deterministic per-protein, so multi-core should
# produce the same ProteinLevelData as single-core: same rows, same
# LogIntensities, same Variance values. The only difference across the two
# paths is how work is distributed to workers.

expect_equal(nrow(QuantDataDefaultLinear$ProteinLevelData),
             nrow(QuantDataParallelLinear$ProteinLevelData),
             info = "Linear multicore should yield same ProteinLevelData row count")

# Sort both by (Protein, RUN) so the comparison is order-independent
linear_single = QuantDataDefaultLinear$ProteinLevelData
linear_multi  = QuantDataParallelLinear$ProteinLevelData
linear_single = linear_single[order(as.character(linear_single$Protein),
                                    as.character(linear_single$RUN)), ]
linear_multi  = linear_multi[order(as.character(linear_multi$Protein),
                                   as.character(linear_multi$RUN)), ]
rownames(linear_single) = NULL
rownames(linear_multi)  = NULL

expect_equal(as.character(linear_single$Protein),
             as.character(linear_multi$Protein),
             info = "Linear multicore should cover the same set of proteins")

expect_equal(linear_single$LogIntensities,
             linear_multi$LogIntensities,
             info = "Linear multicore LogIntensities should match single-core")

if ("Variance" %in% colnames(linear_single) &&
    "Variance" %in% colnames(linear_multi)) {
    expect_equal(linear_single$Variance,
                 linear_multi$Variance,
                 info = "Linear multicore Variance should match single-core")
}


# Test dataProcess with technical replicates & fractions ------------------
msstats_input_fractions_techreps = data.table::fread(
    system.file("tinytest/processed_data/input_techreps_fractions.csv",
                package = "MSstats")
)
QuantDataTechRepsFractions = dataProcess(msstats_input_fractions_techreps,
                                         use_log_file = FALSE)
expect_true(!is.null(QuantDataTechRepsFractions))
expect_true(nrow(QuantDataTechRepsFractions$FeatureLevelData) > 0)

# Test dataProcess with linear summary method and anomaly scores --------------
msstats_input = data.table::data.table(
    ProteinName = c(rep("Q9UFW8", 10), rep("Q96S19", 15)),
    PeptideSequence = c(rep("AEFEEQNVR", 5), rep("TALYVTPLDR", 5),
                        rep("AFPLAEWQPSDVDQR", 5), rep("ASGLLLER", 5),
                        rep("LowAbundancePeptide", 5)),
    PrecursorCharge = rep(2, 25),
    FragmentIon = rep("y3", 25),
    ProductCharge = rep(1, 25),
    IsotopeLabelType = rep("L", 25),
    Condition = rep(c("A", "A", "A", "B", "B"), 5),
    BioReplicate = rep(seq(1:5), 5),
    Run = rep(paste0("Run", seq(1:5)), 5),
    Intensity = c(1000, 1500, 2000, 2500, 3000,
                  1100, 1600, 2100, 2600, 3100,
                  1200, 1700, 2200, 2000, 1700,
                  1300, 1800, 1200, 2800, 1800,
                  100, 200, 300, 400, 500),
    AnomalyScores = c(rep(0.03, 10), c(0.01,0.01,0.01,0.9,0.8), c(0.01,0.01,0.01,0.7,0.5), rep(0.10, 5))
)

# Must run with linear summarization
expect_error(
    dataProcess(msstats_input, summaryMethod="TMP"),
    pattern = "AnomalyScores column detected in your input columns.*set summaryMethod to 'linear'",
    info = "Should error when AnomalyScores is present but summaryMethod is not 'linear'"
)

# Process data with linear summarization
QuantDataLinearAnomaly = dataProcess(msstats_input,
                                     normalization=FALSE, # no normalization
                                     MBimpute=TRUE, # Use imputation
                                     summaryMethod="linear" # Key MSstats+ parameter
)

expect_true("Variance" %in% colnames(QuantDataLinearAnomaly$ProteinLevelData),
            info = "Variance column should be present in ProteinLevelData")
protein_data = QuantDataLinearAnomaly$ProteinLevelData
expect_true(all(protein_data$Variance > 0, na.rm = TRUE),
            info = "All variance values should be positive")
expect_true(all(is.finite(protein_data$Variance), na.rm = TRUE),
            info = "All variance values should be finite")
variance_values = protein_data[protein_data$Protein == "Q9UFW8"]$Variance
expect_true(all(variance_values == variance_values[1]),
            info = paste("Not all variance values are equal for Q9UFW8"))

# Create comparison dataset with low anomaly scores for contrast
msstats_input_low_anomaly = data.table::data.table(
    ProteinName = c(rep("Q9UFW8", 10), rep("Q96S19", 15)),
    PeptideSequence = c(rep("AEFEEQNVR", 5), rep("TALYVTPLDR", 5),
                        rep("AFPLAEWQPSDVDQR", 5), rep("ASGLLLER", 5),
                        rep("LowAbundancePeptide", 5)),
    PrecursorCharge = rep(2, 25),
    FragmentIon = rep("y3", 25),
    ProductCharge = rep(1, 25),
    IsotopeLabelType = rep("L", 25),
    Condition = rep(c("A", "A", "A", "B", "B"), 5),
    BioReplicate = rep(seq(1:5), 5),
    Run = rep(paste0("Run", seq(1:5)), 5),
    Intensity = c(1000, 1500, 2000, 2500, 3000,
                  1100, 1600, 2100, 2600, 3100,
                  1200, 1700, 2200, 2000, 1700,
                  1300, 1800, 1200, 2800, 1800,
                  100, 200, 300, 400, 500),
    AnomalyScores = rep(0.01, 25) # All low anomaly scores
)

QuantDataLowAnomaly = dataProcess(msstats_input_low_anomaly,
                                  normalization=FALSE,
                                  MBimpute=TRUE,
                                  summaryMethod="linear")
high_anomaly_variance = protein_data$Variance
low_anomaly_variance = QuantDataLowAnomaly$ProteinLevelData$Variance
expect_true(mean(high_anomaly_variance, na.rm = TRUE) > mean(low_anomaly_variance, na.rm = TRUE),
            info = "Mean variance should be higher for dataset with high anomaly scores")


# Verify that without AnomalyScores, different behavior occurs
msstats_input_no_anomaly = msstats_input[, !c("AnomalyScores")]
QuantDataNoAnomaly = dataProcess(msstats_input_no_anomaly,
                                 normalization=FALSE,
                                 MBimpute=TRUE,
                                 summaryMethod="linear")
no_anomaly_variance = QuantDataNoAnomaly$ProteinLevelData$Variance
expect_false(identical(high_anomaly_variance, no_anomaly_variance),
             info = "Variance should be different when AnomalyScores are present vs absent")



# Test MSstatsSummarizeSingleLinear function ------------------------------

single_protein = data.table::data.table(
    PROTEIN = c( "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19", "Q96S19" ),
    PEPTIDE = c( "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2", "AFPLAEWQPSDVDQR_2", "ASGLLLER_2", "LowAbundancePeptide_2" ),
    TRANSITION = c( "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1", "y3_1" ),
    FEATURE = c( "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1", "AFPLAEWQPSDVDQR_2_y3_1", "ASGLLLER_2_y3_1", "LowAbundancePeptide_2_y3_1" ),
    LABEL = c( "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L", "L" ),
    GROUP_ORIGINAL = c( "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B" ),
    SUBJECT_ORIGINAL = c( 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5 ),
    RUN = c( "1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5" ),
    GROUP = c( "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B" ),
    SUBJECT = c( "1", "1", "1", "2", "2", "2", "3", "3", "3", "4", "4", "4", "5", "5", "5" ),
    FRACTION = c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ),
    INTENSITY = c( 1200, 1300, 100, 1700, 1800, 200, 2200, 1200, 300, 2000, 2800, 400, 1700, 1800, 500 ),
    ANOMALYSCORES = c( 0.01, 0.01, 0.1, 0.01, 0.01, 0.1, 0.01, 0.01, 0.1, 0.9, 0.7, 0.1, 0.8, 0.5, 0.1 ),
    ABUNDANCE = c( 10.229, 10.344, 6.644, 10.731, 10.814, 7.644, 11.103, 10.229, 8.229, 10.966, 11.451, 8.644, 10.731, 10.814, 8.966 ),
    originalRUN = c( "Run1", "Run1", "Run1", "Run2", "Run2", "Run2", "Run3", "Run3", "Run3", "Run4", "Run4", "Run4", "Run5", "Run5", "Run5" ),
    censored = c( FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE ),
    newABUNDANCE = c( 10.229, 10.344, NA, 10.731, 10.814, NA, 11.103, 10.229, NA, 10.966, 11.451, NA, 10.731, 10.814, NA ),
    cen = c( 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0 ),
    nonmissing = c( TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE ),
    n_obs = c( 5, 5, 0, 5, 5, 0, 5, 5, 0, 5, 5, 0, 5, 5, 0 ),
    n_obs_run = c( 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 ),
    total_features = c( 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 ),
    prop_features = c( 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667, 0.667 ),
    nonmissing_all = c( TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE ),
    ABUNDANCE_cut = c( 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA, 10.127, 10.127, NA ),
    any_censored = c( FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE )
)

summarized = MSstats::MSstatsSummarizeSingleLinear(
    single_protein,
    impute = FALSE,
    censored_symbol = "NA",
    remove50missing = FALSE,
    aft_iterations = 90,
    equal_variances = TRUE
)

protein_level_summary = summarized[[1]]
feature_level_summary = summarized[[2]]

## Basic tests
expect_equal(nrow(protein_level_summary), 5)
expect_equal(ncol(protein_level_summary), 4)
expect_true(all(protein_level_summary$Protein == "Q96S19"))
expect_equal(as.numeric(protein_level_summary$RUN), 1:5)
expect_equal(nrow(feature_level_summary), 10)
expect_equal(ncol(feature_level_summary), 5)
expect_equal(length(unique(feature_level_summary$FEATURE)), 2)
expect_equal(length(unique(feature_level_summary$RUN)), 5)

## Verify that variance is higher in runs 4-5 with high anomaly scores
variance_runs_1_3 = protein_level_summary[RUN %in% 1:3, Variance]
variance_runs_4_5 = protein_level_summary[RUN %in% 4:5, Variance]

for (var_45 in variance_runs_4_5) {
    for (var_13 in variance_runs_1_3) {
        expect_true(var_45 > var_13,
                  info = paste("Variance in runs 4-5 should be > variance in runs 1-3:",
                               var_45, ">", var_13))
    }
}

mean_variance_1_3 = mean(variance_runs_1_3)
mean_variance_4_5 = mean(variance_runs_4_5)

expect_true(mean_variance_4_5 > 10 * mean_variance_1_3,
          info = paste("Mean variance of runs 4-5 should be at least 10x higher:",
                       mean_variance_4_5, "vs", 10 * mean_variance_1_3))

## Verify that protein levels in runs 4-5 are closer to ASGLLLER than AFPLAEWQPSDVDQR
protein_runs_4_5 = protein_level_summary[RUN %in% 4:5]
afpla_runs_4_5 = feature_level_summary[FEATURE == "AFPLAEWQPSDVDQR_2_y3_1" & RUN %in% 4:5]
asgll_runs_4_5 = feature_level_summary[FEATURE == "ASGLLLER_2_y3_1" & RUN %in% 4:5]
for (i in 1:nrow(protein_runs_4_5)) {
    run_num = protein_runs_4_5$RUN[i]
    protein_logint = protein_runs_4_5$LogIntensities[i]

    afpla_abundance = afpla_runs_4_5[RUN == run_num, newABUNDANCE]
    asgll_abundance = asgll_runs_4_5[RUN == run_num, newABUNDANCE]

    dist_to_afpla = abs(protein_logint - afpla_abundance)
    dist_to_asgll = abs(protein_logint - asgll_abundance)

    expect_true(dist_to_asgll < dist_to_afpla,
              info = paste("Run", run_num, ": Protein LogIntensity", protein_logint,
                           "should be closer to ASGLLLER", asgll_abundance,
                           "(dist =", round(dist_to_asgll, 4), ")",
                           "than to AFPLAEWQPSDVDQR", afpla_abundance,
                           "(dist =", round(dist_to_afpla, 4), ")"))
}


# =============================================================================
# Integration tests: dataProcess memory efficiency
# =============================================================================
#
# One in-process measurement captures both numbers we care about per run:
#
#   - retained_mb: Vcells "used (Mb)" delta after the call + gc(). Catches
#     output-size bloat and objects that remain reachable after return
#     (leaked references in package options, closures, or the global env).
#
#   - peak_mb: Vcells "max used (Mb)" delta across the call, read from gc()
#     after gc(reset = TRUE). R updates max-used on every allocation, so
#     transient spikes that are freed before return — merge copies,
#     as.data.frame deep copies, per-protein lm/survreg model objects —
#     still show up here.
#
# Scope limit: gc() only sees the current R process. In the parallel
# (numberOfCores > 1) path, worker processes run in their own address spaces
# and are invisible to this measurement. peak_mb there reflects the parent
# only (serialization buffers, result collation); it will be much smaller
# than total system memory pressure, but still catches any parent-side
# regression. An earlier callr/ps-based subprocess sampler tried to see
# workers, but dataProcess under pkgload::load_all in a fresh subprocess
# did not finish in bounded time, so that approach was removed.
# =============================================================================

# --- Helper: measure retained + peak memory for a single call ---------------
#
# Returns a list with:
#   $result          return value of expr_fn
#   $input_mb        size of input_data in MB (for ratio reporting)
#   $retained_mb     Vcells used delta after expr_fn + gc()
#   $peak_mb         Vcells max-used delta seen during expr_fn
#   $retained_ratio  retained_mb / input_mb
#   $peak_ratio      peak_mb    / input_mb
#   $elapsed_s       wall-clock time of expr_fn
#
# gc()'s matrix layout varies by R version (some versions add a "limit (Mb)"
# column), so we find the "(Mb)" column that follows each labelled count
# column by name lookup rather than hardcoding positions.
measure_memory = function(expr_fn, input_data) {
    input_mb = as.numeric(object.size(input_data)) / 1e6
    gc()                       # free any dead objects so baseline is clean
    gc(reset = TRUE)           # reset max-used tracking to current used
    before = gc()
    cols = colnames(before)
    used_mb_col = which(cols == "used")[1] + 1      # (Mb) right after "used"
    max_mb_col  = which(cols == "max used")[1] + 1  # (Mb) right after "max used"
    baseline_mb = before["Vcells", used_mb_col]
    t0 = Sys.time()
    result = expr_fn()
    elapsed_s = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    gc()                       # free intermediates for retained measurement
    after = gc()               # post-gc state also carries max-used
    retained_mb = after["Vcells", used_mb_col] - baseline_mb
    peak_mb     = after["Vcells", max_mb_col]  - baseline_mb
    list(result = result,
         input_mb = input_mb,
         retained_mb = retained_mb,
         peak_mb = peak_mb,
         retained_ratio = retained_mb / input_mb,
         peak_ratio = peak_mb / input_mb,
         elapsed_s = elapsed_s)
}


# --- Build a test dataset large enough for meaningful memory measurement ----
#
# SRMRawData (720 rows) and DDARawData (2070 rows) are too small to see
# memory differences clearly against GC noise. Replicate DDARawData 100x,
# synthesizing unique ProteinName suffixes so dataProcess treats each
# replicate as a separate protein.
set.seed(42)
base_data = data.table::as.data.table(DDARawData)
n_replicates = 100
replicated_data = data.table::rbindlist(lapply(seq_len(n_replicates), function(i) {
    d = data.table::copy(base_data)
    d$ProteinName = paste0(d$ProteinName, "_rep", i)
    d
}))
replicated_data = as.data.frame(replicated_data)
input_size_mb = as.numeric(object.size(replicated_data)) / 1e6

cat(sprintf("Memory test dataset: %d rows, %.1f MB\n",
            nrow(replicated_data), input_size_mb))


# Shared helper for reporting + asserting on a measure_memory() result.
# retained_limit and peak_limit are ratios vs input_mb.
check_memory = function(label, mem, retained_limit, peak_limit) {
    cat(sprintf(
        "%s: input=%.1f MB, retained=%.1f MB (%.1fx), peak=%.1f MB (%.1fx), %.1fs\n",
        label, mem$input_mb,
        mem$retained_mb, mem$retained_ratio,
        mem$peak_mb,     mem$peak_ratio,
        mem$elapsed_s))
    expect_true(mem$retained_ratio < retained_limit,
                info = sprintf(
                    "%s retained memory should be < %gx input. retained=%.1f MB, ratio=%.1fx",
                    label, retained_limit, mem$retained_mb, mem$retained_ratio))
    expect_true(mem$peak_ratio < peak_limit,
                info = sprintf(
                    "%s peak memory should be < %gx input. peak=%.1f MB, ratio=%.1fx",
                    label, peak_limit, mem$peak_mb, mem$peak_ratio))
    expect_true(nrow(mem$result$FeatureLevelData) > 0,
                info = sprintf("%s: FeatureLevelData should have rows", label))
    expect_true(nrow(mem$result$ProteinLevelData) > 0,
                info = sprintf("%s: ProteinLevelData should have rows", label))
}


MSstatsConvert::MSstatsLogsSettings(FALSE)

# Thresholds below are calibrated against the current post-fix baseline on
# this specific 100-replicate fixture (~10.8 MB input). Peak Vcells sits
# around 55-60x input here because the fixture is tiny — fixed overhead
# (factor levels, model infrastructure) dominates. The retained ratio is a
# more portable regression signal and stays around 1.5x.
#
# Peak memory grows SUPERLINEARLY with input on the current implementation
# (measured: n=100 -> 615 MB, n=300 -> 4.6 GB, n=500 -> 12.4 GB on this
# machine). This test only guards the small-fixture behaviour; it is NOT a
# statement that peak memory scales safely to production-sized data. Tuning
# the pipeline so peak stays sub-linear is tracked separately.
RETAINED_LIMIT = 5    # output is ~1x input; 5x catches leaked references
PEAK_LIMIT_TMP = 75   # post-fix baseline on 10 MB fixture ~57x; 75x leaves
                      # ~30% slack for R/OS variance without going flaky
PEAK_LIMIT_LINEAR = 75

# Single-core TMP (default).
single_core = measure_memory(
    function() dataProcess(replicated_data, use_log_file = FALSE),
    replicated_data
)
check_memory("Single-core TMP", single_core,
             retained_limit = RETAINED_LIMIT, peak_limit = PEAK_LIMIT_TMP)
rm(single_core); gc()

# Multi-core TMP (2 workers).
# peak_mb here is parent-process only — worker processes are invisible to
# gc(). This is a regression guard on parent memory only; serialization
# buffers + result collation should stay within the same budget as the
# single-core parent.
multi_core = measure_memory(
    function() dataProcess(replicated_data, use_log_file = FALSE,
                           numberOfCores = 2),
    replicated_data
)
check_memory("Multi-core TMP (2)", multi_core,
             retained_limit = RETAINED_LIMIT, peak_limit = PEAK_LIMIT_TMP)
rm(multi_core); gc()

# Linear single-core.
# Linear path uses lm + survreg per protein; watches for Issue 1 regressions
# where model objects retain qr/residuals/model frame past their last use.
linear_mem = measure_memory(
    function() dataProcess(replicated_data, use_log_file = FALSE,
                           summaryMethod = "linear"),
    replicated_data
)
check_memory("Linear single-core", linear_mem,
             retained_limit = RETAINED_LIMIT, peak_limit = PEAK_LIMIT_LINEAR)
rm(linear_mem); gc()

# Linear multi-core (2 workers).
# Peak_mb here is parent-process only — workers are invisible to gc() — so
# this is a regression guard on parent memory only. Pairs with the
# "Linear single-core" test above: both should land at the same budget,
# because the parent's job in the parallel path (serialize per-protein
# chunks, collate results) is the same size as the single-core loop.
linear_multi_mem = measure_memory(
    function() dataProcess(replicated_data, use_log_file = FALSE,
                           summaryMethod = "linear", numberOfCores = 2),
    replicated_data
)
check_memory("Linear multi-core (2)", linear_multi_mem,
             retained_limit = RETAINED_LIMIT, peak_limit = PEAK_LIMIT_LINEAR)
rm(linear_multi_mem); gc()


# =============================================================================
# profmem diagnostic — no single allocation should exceed 150 MB
# =============================================================================
#
# Catches a different failure than peak: a regression that allocates one
# large object (e.g. the pre-fix .finalizeTMP merge at ~350 MB) is pinned
# directly to its allocation call, making it trivially diagnosable. Skipped
# cleanly if profmem is unavailable. Single-core only — Rprofmem doesn't see
# forked worker processes.

if (!requireNamespace("profmem", quietly = TRUE)) {
    exit_file("profmem not installed — skipping allocation diagnostic")
}

pm = profmem::profmem(
    dataProcess(replicated_data, use_log_file = FALSE)
)
big_allocs = pm[!is.na(pm$bytes) & pm$bytes > 150 * 1024^2, ]
expect_equal(nrow(big_allocs), 0L,
             info = paste("Single allocation(s) >150 MB detected:",
                          paste(utils::capture.output(print(big_allocs)),
                                collapse = "\n")))
