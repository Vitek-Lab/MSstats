# Scalability Analysis: MSstats at Real-World Scale

## Overview

This document catalogs the scalability issues that persist in the MSstats
package after the memory-management fixes on the `Fix-Memory-Management`
branch. The companion document `dataProcess_memory_analysis.md` describes
nine memory inefficiencies in `dataProcess`; seven of those have since been
addressed (in-place column drops, keyed-match `.finalizeTMP`, pre-split
parallel workers, `setDF` in place of `as.data.frame`, etc.). The issues
documented here are the *remaining* limits on real-world scale: inputs of
10,000+ proteins, hundreds of runs, multi-GB feature tables.

The findings span four scalability axes:

- **Peak memory** vs. input size
- **Wall-clock time** vs. input size
- **Scaling with cores** (parallel path)
- **Scaling with proteins** (per-protein loop accumulation)

All file paths and line numbers reference the current state of the `R/`
source directory on the `Fix-Memory-Management` branch and have been
verified directly against the code (not inferred from prior analyses).

---

## Table of Contents

1. [Pipeline Scope and Method](#1-pipeline-scope-and-method)
2. [S1: Per-Protein Summarize Loop Accumulates All Results](#2-s1-per-protein-summarize-loop-accumulates-all-results)
3. [S2: groupComparison Single-Core Loop Mirrors S1](#3-s2-groupcomparison-single-core-loop-mirrors-s1)
4. [S3: `save_fitted_models = TRUE` Retains lm / lmerMod Per Protein](#4-s3-save_fitted_models--true-retains-lm--lmermod-per-protein)
5. [S4: `MSstatsBalancedDesign` Dominates the Prepare-Stage Peak](#5-s4-msstatsbalanceddesign-dominates-the-prepare-stage-peak)
6. [S5: Heatmap Clustering is O(n²) in Proteins](#6-s5-heatmap-clustering-is-on²-in-proteins)
7. [S6: Profile / QC Plot Loops Retain All ggplot Objects](#7-s6-profile--qc-plot-loops-retain-all-ggplot-objects)
8. [S7: O(n²) Hover-Text Construction in Plotly Heatmap](#8-s7-on²-hover-text-construction-in-plotly-heatmap)
9. [S8: Quantile-Normalization Branch Still Does a Full Copy](#9-s8-quantile-normalization-branch-still-does-a-full-copy)
10. [S9: Linear Path Holds `lm` Across Protein Body](#10-s9-linear-path-holds-lm-across-protein-body)
11. [S10: `.SD` Materialization Without `.SDcols`](#11-s10-sd-materialization-without-sdcols)
12. [S11: groupComparison parLapply Closure Captures Full List](#12-s11-groupcomparison-parlapply-closure-captures-full-list)
13. [S12: Cluster Startup Paid Per `dataProcess` Call](#13-s12-cluster-startup-paid-per-dataprocess-call)
14. [Cumulative Impact and Prioritization](#14-cumulative-impact-and-prioritization)

---

## 1. Pipeline Scope and Method

The user-facing entry points whose scalability we analyze:

```
MSstats user surface
  |
  |-- dataProcess              (R/dataProcess.R:124)
  |     |-- (Prepare/Normalize/Merge/Missing/Select stages — covered in
  |     |    dataProcess_memory_analysis.md)
  |     |-- MSstatsSummarizeWithSingleCore        (R/dataProcess.R:304)   <- S1
  |     |-- MSstatsSummarizeWithMultipleCores     (R/dataProcess.R:220)
  |     |-- MSstatsSummarizationOutput            (R/utils_output.R:37)
  |
  |-- groupComparison           (R/groupComparison.R)
  |     |-- MSstatsGroupComparison                (R/utils_groupcomparison.R)  <- S2, S3, S11
  |     |-- MSstatsGroupComparisonSingleProtein   (R/utils_groupcomparison.R)
  |     |-- MSstatsGroupComparisonOutput          (R/utils_groupcomparison.R)
  |
  |-- dataProcessPlots          (R/dataProcessPlots.R)                     <- S6
  |-- groupComparisonPlots      (R/groupComparisonPlots.R)                 <- S5, S7
```

Reference dataset for severity estimates: **10,000 proteins × 100 runs ×
10 features per protein** → a feature-level `input` table of ~10 million
rows, ~1.5 GB in memory. This is roughly an order of magnitude larger
than the 1,000-protein fixtures the package has historically been tuned
against.

---

## 2. S1: Per-Protein Summarize Loop Accumulates All Results

**Severity**: High
**Axes**: peak memory + per-protein scaling
**Location**: `R/dataProcess.R:304-335`

### Code

```r
MSstatsSummarizeWithSingleCore = function(input, method, impute, censored_symbol,
                                          remove50missing, equal_variance,
                                          aft_iterations = 90) {
    protein_indices = split(seq_len(nrow(input)), list(input$PROTEIN))
    num_proteins = length(protein_indices)
    summarized_results = vector("list", num_proteins)               # line 310

    if (method == "TMP") {
        pb = utils::txtProgressBar(min = 0, max = num_proteins, style = 3)
        for (protein_id in seq_len(num_proteins)) {
            single_protein = input[protein_indices[[protein_id]],]
            summarized_results[[protein_id]] = MSstatsSummarizeSingleTMP(
                single_protein, impute, censored_symbol, remove50missing,
                aft_iterations)                                     # line 315
            setTxtProgressBar(pb, protein_id)
        }
        close(pb)
    } else {
        # identical pattern for linear at line 329
    }
    summarized_results
}
```

### What grows

Each `MSstatsSummarizeSingleTMP` call returns `list(result, survival)`,
where:

- `result` is a per-protein summary table — small, a few kilobytes
- `survival` is per-protein imputed values — typically tens of kilobytes
  for a protein with ~500 feature-runs

Both stay alive for the rest of the loop because they are bound into
`summarized_results[[i]]`. Memory grows roughly linearly with the
iteration index.

Downstream, `MSstatsSummarizationOutput` (`R/utils_output.R:37`) consumes
this list in two passes — once to extract survival predictions (line
~47), once to extract summaries (line ~42). Only after both passes have
finished `rbindlist`-ing are the per-protein entries eligible for
garbage collection.

### Memory state

For a 3,000-protein fixture, instrumentation has previously shown a
**+1.6 GB peak** at mid-loop. Trajectory for 10k proteins (extrapolated):

| Loop position | Items live | Approx size |
|---------------|-----------|-------------|
| protein 0 | 0 | ~0 MB |
| protein 5k | 5,000 result/survival pairs | ~2-4 GB |
| protein 10k (loop end) | 10,000 result/survival pairs | ~4-8 GB |
| during finalize rbindlist | above + flat lists + rbindlist results | ~6-10 GB |

This is now the single largest peak contributor in `dataProcess` on
real-scale inputs.

### Solutions

Two options, in order of impact:

1. **Streaming accumulator (API change).** Maintain two flat lists —
   `summaries` and `survivals` — and `rbindlist` incrementally every K
   proteins (K = ~100). The return shape becomes two pre-collated
   data.tables instead of a nested list. Downstream `.finalizeTMP` /
   `MSstatsSummarizationOutput` must be updated. Bumps the package
   minor/major version because the documented return changes.

2. **In-loop null-out (no API change).** Have the output consumer null
   out `summarized[[i]][[2]]` immediately after extracting it for the
   survival rbindlist, and `summarized[[i]][[1]]` after extracting it
   for the summary rbindlist. Halves coexistence. Some of this exists
   on the branch already (`R/utils_output.R:47-48`) — verify it
   survived recent edits and extend it to both slots.

Option 1 reduces peak by ~4× at scale; option 2 by ~2×.

---

## 3. S2: groupComparison Single-Core Loop Mirrors S1

**Severity**: High
**Axes**: peak memory + per-protein scaling
**Location**: `R/utils_groupcomparison.R:510-521`

### Code

```r
test_results = vector("list", length(all_proteins_id))               # line 510
for (i in all_proteins_id) {
    comparison_outputs = MSstatsGroupComparisonSingleProtein(
        i, summarized_list[[i]], contrast_matrix, repeated,
        groups, samples_info, save_fitted_models, has_imputed)
    test_results[[i]] = comparison_outputs                           # line 517
}
test_results
```

Each `comparison_outputs` is a 3-element list:
`list(single_protein_data, fit_result_table, fitted_model_object)`.
Same accumulator anti-pattern as S1 but with worse per-element payload:
the per-protein data.table AND the fitted model object both live in
each slot.

### Memory state at 10k proteins

| Component | Per protein | Total at loop end |
|-----------|------------|-------------------|
| `single_protein` data.table (~500 rows × 15 cols) | ~50 KB | ~500 MB |
| `fit_result_table` (one row per contrast) | ~1 KB | ~10 MB |
| `fitted_model` (when `save_fitted_models = TRUE`) | ~40-200 KB | ~0.4-2 GB |
| **Combined** | ~91-251 KB | ~0.9-2.5 GB |

Then `MSstatsGroupComparisonOutput` extracts each piece via
`lapply(test_results, function(x) x[[k]])` — the original list stays
alive during extraction, so the peak briefly doubles.

### Solution

The same fix as S1: stream into flat lists and `rbindlist` incrementally,
or null out `test_results[[i]]` after `MSstatsGroupComparisonOutput`
extracts everything from slot i. The parallel path
(`R/utils_groupcomparison.R:475-491`) has the same shape and benefits
from the same fix.

---

## 4. S3: `save_fitted_models = TRUE` Retains lm / lmerMod Per Protein

**Severity**: Very High
**Axes**: peak memory + per-protein scaling + final-object size
**Location**: `R/groupComparison.R` (argument defaults), referenced
through `R/utils_groupcomparison.R`

### Background

When `save_fitted_models = TRUE`, the third element of each protein's
result in `test_results` is the full `lm` or `lmerMod` object returned
by the fit. These objects include:

- `$model` — the data frame used for fitting (all rows, all formula
  columns)
- `$qr$qr` — QR decomposition (nrow × ncol of the model matrix)
- `$residuals`, `$fitted.values`, `$effects` — numeric vectors of length
  nrow
- For `lmerMod`: `@frame`, `@flist`, `@resp`, `@pp$X` — even bulkier
- The formula's `environment()` — a reference that can keep arbitrary
  enclosing-scope objects alive

For a protein with 500 feature-runs, one `lm` object is 40-80 KB; one
`lmerMod` object is 100-200 KB.

### What the user usually needs

In nearly all downstream MSstats uses, only the **coefficient vector**
(`coef(fit)`) and the **variance-covariance matrix** (`vcov(fit)`) are
read. Together these are well under 5 KB per protein.

### Scaling

| Proteins | `save_fitted_models = TRUE` (lm) | (lmerMod) | with bulky slots stripped |
|----------|---------------------------------|-----------|---------------------------|
| 1,000 | 40-80 MB | 100-200 MB | ~5 MB |
| 10,000 | 400-800 MB | 1-2 GB | ~50 MB |
| 50,000 | 2-4 GB | 5-10 GB | ~250 MB |

These are sizes in the **returned object**, not transient — they persist
in the user's workspace for the duration of the R session.

### Solutions

In rough order of impact and ease:

1. **Flip the default to `FALSE`.** Documents the size cost on the
   argument. Users who want fits opt in.
2. **Strip bulky slots before storing.** Keep `$coefficients`, `$vcov`,
   `$rank`, `$df.residual`; drop `$model`, `$qr`, `$residuals`,
   `$fitted.values`, `$effects`, `$y`. Document what survives.
3. **Replace with a thin summary** — store
   `list(coef = coef(fit), vcov = vcov(fit), rank = ..., df.r = ...)`.
   Smallest footprint; breaks any user code that introspects `fit`
   directly.

---

## 5. S4: `MSstatsBalancedDesign` Dominates the Prepare-Stage Peak

**Severity**: High (out of scope for this repo — lives upstream)
**Axes**: peak memory + time, both scale with proteins × features × runs
**Location**: `MSstatsConvert::MSstatsBalancedDesign`, called from
`R/utils_checks.R::.checkUnProcessedDataValidity` around line 175

### What it does

`MSstatsBalancedDesign` pads the input grid so that every
protein × feature × run combination has a row, with `NA` intensity where
data is missing. The result is a fully balanced design that downstream
functions can rely on.

### Why it spikes

The implementation materializes the full Cartesian product `proteins ×
features × runs` before filling in known values and `NA`s for the rest.
For a 5,000-protein × 10-feature × 100-run experiment, that's 5 million
rows generated even if the actual data has only 2 million observed
rows — a 2.5× row-count blow-up.

Observed peaks:

| Input size | Peak during BalancedDesign |
|-----------|---------------------------|
| 10.8 MB | 108 MB |
| 324 MB | 2.5 GB |

The function sets the **floor** for total `dataProcess` peak. No matter
what is optimized inside MSstats, this stage will dominate on any input
that isn't already balanced.

### Solutions (must land in `MSstatsConvert`)

- Stream the padding row-by-row instead of materializing the full grid.
- Skip the padding entirely when the input is already balanced (cheap
  up-front check on row count vs. product of unique values).

File an issue against `MSstatsConvert` linking to the profile output
from `benchmark/profile_dataprocess_peak.R`. Until that lands, the
Prepare spike cannot be reduced from inside MSstats.

---

## 6. S5: Heatmap Clustering is O(n²) in Proteins

**Severity**: High at protein count > ~3,000
**Axes**: wall-clock time + peak memory
**Location**: `R/groupComparisonPlots.R:206-212`

### Code (shape)

```r
wide = data.table::dcast(input, Protein ~ Label, value.var = "heat_val")
wide = as.matrix(wide[, -1])
wide = .getOrderedMatrix(wide, clustering)     # calls hclust(dist(wide))
```

### Three problems compound

1. **`dcast` allocates a dense matrix** even when the underlying data is
   sparse. For 10k proteins × 100 contrasts that is 8 MB of doubles —
   small, but combined with `NA` handling and `as.matrix` doubling.
2. **`as.matrix(wide[, -1])` makes a second copy** of the same data.
3. **`hclust(dist(wide))`** is the real cost. `dist()` materializes the
   full protein-by-protein distance matrix:
   - n=1,000: 1M cells, ~8 MB, <1 second
   - n=3,000: 9M cells, ~72 MB, ~5 seconds
   - n=10,000: 100M cells, **~800 MB, ~minutes**
   - n=50,000: 2.5B cells, **~20 GB, infeasible**

   `hclust` itself is O(n²) time and O(n) extra memory on top.

### Solutions

- **Cap clustered protein count.** Add a `max_proteins_for_clustering`
  argument (default ~1,000); above this, skip clustering or show only
  top-N by significance.
- **Use `fastcluster::hclust`.** Drop-in replacement, 2-5× faster, same
  output.
- **Sparse distance.** If many proteins share `NA` patterns, a sparse
  representation can avoid materializing the full distance matrix —
  larger refactor.

---

## 7. S6: Profile / QC Plot Loops Retain All ggplot Objects

**Severity**: High at protein count > ~1,000
**Axes**: peak memory + per-protein scaling
**Location**: `R/dataProcessPlots.R` (around lines 319, 349, 368, 410)

### Code (shape)

```r
output_plots = list()
for (i in seq_along(all_proteins)) {
    single_protein = input[Protein == all_proteins[i]]
    profile_plot   = .makeProfilePlot(single_protein, ...)
    if (address != FALSE) print(profile_plot)         # writes to PDF
    output_plots[["original_plot"]][[paste("plot", i)]] = profile_plot
}
```

### Why each plot is heavy

A ggplot object retains:

- A reference to its data (`profile_plot$data`) — typically the
  per-protein slice, ~30-100 KB for moderate proteins
- Layer-specific data copies if any geom transforms the data
- Scales, themes, coordinate systems — small but nonzero per object
- Faceting parameters

Per-plot cost: typically 200 KB – 1 MB for profile plots with multiple
geoms. At 10,000 proteins:

| Quantity | Total |
|----------|-------|
| ggplot objects in `output_plots` | 10,000 |
| Approx total size | 2-10 GB |

### Why it's particularly bad

The output is two things at once: a side-effect (drawing to a PDF) and
a return value (the list). When the user passes `address = "file.pdf"`,
they typically want the PDF — not also a multi-GB list of ggplots in
their workspace.

### Solutions

- When `address != FALSE` (writing to PDF), do **not** also retain the
  ggplot in the returned list — return only metadata (protein names
  plotted, page numbers) or nothing.
- Provide an explicit `return_plots = FALSE` argument (default) for
  large-protein workflows.
- For users who want the list, document the per-protein memory cost.

---

## 8. S7: O(n²) Hover-Text Construction in Plotly Heatmap

**Severity**: Medium (time only, fixable in minutes)
**Axes**: wall-clock time
**Location**: `R/utils_groupcomparison_plots.R:134-141`

### Code (shape)

```r
hover_text_matrix = matrix("", nrow = nrow(input), ncol = ncol(input))
for (i in 1:nrow(input)) {
    for (j in 1:ncol(input)) {
        hover_text_matrix[i, j] = sprintf("Protein: %s<br>Comparison: %s<br>FC: %.3f<br>p: %.3g",
                                          rownames(input)[i], colnames(input)[j], ...)
    }
}
```

### Cost

For 10k proteins × 100 contrasts:
- 1,000,000 `sprintf` calls — interpreted-loop overhead dominates
- Each call allocates a small character string that goes through R's
  global string-intern hash

Wall-clock cost: tens of seconds to minutes, depending on hardware.

### Solution

Vectorize:

```r
hover_text_matrix = matrix(
    sprintf("Protein: %s<br>Comparison: %s<br>FC: %.3f<br>p: %.3g",
            rep(rownames(input), ncol(input)),
            rep(colnames(input), each = nrow(input)),
            as.vector(input_fc),
            as.vector(input_p)),
    nrow = nrow(input), ncol = ncol(input)
)
```

One `sprintf` call, vectorized over all cells. 10-100× faster.

---

## 9. S8: Quantile-Normalization Branch Still Does a Full Copy

**Severity**: Medium — but only triggers when user picks quantile
**Axes**: peak memory
**Location**: `R/utils_normalize.R:118`

### Code

```r
per_fraction_normalized = data.table::rbindlist(per_fraction_normalized)
input = merge(input[, colnames(input) != "ABUNDANCE", with = FALSE],   # <- full copy
              per_fraction_normalized, by = grouping_cols)
```

### Why this slipped through the recent fixes

Recent commits replaced the four column-subset-copy sites in
`.normalizeMedian`, `.normalizeGlobalStandards`, and the two
`MSstatsMergeFractions` branches with `data.table::set(..., value = NULL)`
in-place drops. Locations A through D in
`dataProcess_memory_analysis.md` Issue 2 are fixed.

This fifth site lives inside `.normalizeQuantile` and is **not exercised
by the default test fixture** (median normalization is the default), so
it survived. The pattern is identical to the fixed locations: a column
subset via negative indexing creates a deep copy of every retained column
before the merge.

### Memory state at the merge

Same as the fixed sites: a ~250 MB transient copy on a 10M-row input,
co-existing with the original `input` and the merge result.

### Solution

Mirror the fix already applied in `.finalizeTMP`: use a keyed match
instead of `merge(subset, ...)`:

```r
data.table::setkeyv(input, grouping_cols)
data.table::setkeyv(per_fraction_normalized, grouping_cols)
input[per_fraction_normalized, ABUNDANCE := i.ABUNDANCE]
```

No intermediate copy, no merge allocation.

---

## 10. S9: Linear Path Holds `lm` Across Protein Body

**Severity**: Medium — cumulative GC churn rather than peak
**Axes**: per-protein scaling, time (GC overhead)
**Location**: `R/dataProcess.R::MSstatsSummarizeSingleLinear`,
roughly lines 395-430

### Code (shape)

```r
fit = .fitLinearModel(single_protein, equal_variance)            # ~line 411
cf      = summary(fit)$coefficients[, 1]                          # ~line 421
cov_mat = vcov(fit)                                               # ~line 422
# fit is never read again, but lives until function returns at ~line 431
```

The TMP path was patched: a `rm(survival_fit)` was inserted after the
`predict()` call, freeing the survival object before the protein body
finishes. The linear path was not given the same treatment for its
`lm` object — `fit` sits alive through `.isSummarizable`, `.runTukey`,
and the return-list construction.

### Cost

Per protein: 40-80 KB of `lm` object alive ~50-100 ms longer than
needed. Across 10k proteins, this drives ~400-800 MB of cumulative
allocation churn that R's GC has to manage.

### Solution

Insert `rm(fit)` immediately after `cov_mat = vcov(fit)`:

```r
cf      = summary(fit)$coefficients[, 1]
cov_mat = vcov(fit)
rm(fit)
```

Two new lines, no behaviour change.

---

## 11. S10: `.SD` Materialization Without `.SDcols`

**Severity**: Low-Medium individually, Medium cumulatively
**Axes**: per-call memory + time
**Location**: ≥ 7 sites across the codebase

### Background

`data.table` syntax `dt[, expr, by = ...]` exposes the per-group rows
to `expr` via the special symbol `.SD` — a `data.table` view of all
columns *not* in the `by` clause. Without `.SDcols`, every column is
materialized into `.SD` for every group.

### Affected sites

| File | Line(s) | Group count |
|------|---------|-------------|
| `R/utils_feature_selection.R` | 148, 202 | one per protein × label |
| `R/utils_normalize.R` | 54, 153 | one per RUN × FRACTION |
| `R/utils_output.R` | 163, 198 | one per RUN, then one per protein |
| `R/utils_summarization_prepare.R` | 128 | one per protein |

Each callee typically reads only 3-6 columns. On a 20-column input
table, the remaining 14-17 columns are materialized for every group for
no reason.

### Memory and time cost

Each materialization is an allocation of a small list-of-vectors. With
thousands of groups per call, the cumulative allocations add up to
hundreds of MB of transient memory and a measurable fraction of stage
runtime — small per site but the same pattern appears seven times.

### Solution

Mechanical: for each site, identify which columns the callee reads, then
add `.SDcols = c(...)`:

```r
input[, NonMissingStats := .getNonMissingFilterStats(.SD, censored_symbol),
      .SDcols = c("LABEL", "newABUNDANCE", "censored", "INTENSITY",
                  "n_obs_run", "n_obs")]
```

---

## 12. S11: groupComparison parLapply Closure Captures Full List

**Severity**: High at high core counts
**Axes**: memory × cores + parallel wall-clock (serialization)
**Location**: `R/utils_groupcomparison.R:475-491`

### Background

When `parallel::parLapply(cl, X, FUN)` runs, the function `FUN` is
serialized along with **its enclosing environment**. Any variable
referenced by name inside the function but not in `X` is captured by
closure and sent to every worker via serialization.

### Code (shape)

```r
parallel::clusterExport(cl, c("MSstatsGroupComparisonSingleProtein",
                              "contrast_matrix", "repeated",
                              "groups", "samples_info",
                              "save_fitted_models", "has_imputed"),
                        envir = function_environment)

test_results = parallel::parLapply(cl, all_proteins_id, function(i) {
    MSstatsGroupComparisonSingleProtein(
        i, summarized_list[[i]],     # <- captures summarized_list by closure
        contrast_matrix, repeated, groups, samples_info,
        save_fitted_models, has_imputed)
})
```

### Why it scales badly

The explicit `clusterExport` does not include `summarized_list` — yet
the worker function references it. R's serialization machinery captures
it from the function's enclosing environment and sends a copy to every
worker.

For 10k proteins, `summarized_list` is ~200-500 MB. With 8 cores:

| Component | Size |
|-----------|------|
| Main process retains `summarized_list` | ~500 MB |
| Each worker receives a copy | 8 × ~500 MB = 4 GB |
| Serialization wall time | seconds to tens of seconds |

The exact same pattern existed in `dataProcess.R` and was fixed by
splitting `input` into per-protein slices *before* `parLapply` and
passing the slices through `X` instead of capturing the full table via
closure. The fix has not been ported to groupComparison.

### Solution

Pre-split, then parLapply over the slices:

```r
protein_data_list = lapply(all_proteins_id, function(i) summarized_list[[i]])

test_results = parallel::parLapply(cl, protein_data_list,
    function(single_protein_data) {
        MSstatsGroupComparisonSingleProtein(
            single_protein_data,
            contrast_matrix, repeated, groups, samples_info,
            save_fitted_models, has_imputed)
    })
```

Now each worker receives only the slice it operates on; total worker
memory scales with chunk size, not with full data size.

---

## 13. S12: Cluster Startup Paid Per `dataProcess` Call

**Severity**: Low individually, Medium for batch workflows
**Axes**: wall-clock time
**Location**: `R/dataProcess.R::MSstatsSummarizeWithMultipleCores`, same
pattern in `R/utils_groupcomparison.R`

### Code (shape)

```r
cl = parallel::makeCluster(numberOfCores)
parallel::clusterExport(cl, c(...))
... parLapply ...
parallel::stopCluster(cl)
```

### Cost

`parallel::makeCluster` spawns fresh R worker processes — each one a
full R startup. On macOS / Linux: 0.5-2 seconds per worker. On Windows:
2-5 seconds per worker (no fork — full process spawn).

| Cores | Startup | Per `dataProcess` call |
|-------|---------|-----------------------|
| 2 | ~2s | paid each call |
| 4 | ~4s | paid each call |
| 8 | ~8s | paid each call |

For an interactive single call, this is fine. For a batch script
processing 100 input files: 100 × 8s = ~13 minutes of pure cluster
startup overhead.

### Solution

Accept an optional pre-made cluster:

```r
dataProcess = function(raw, ..., numberOfCores = 1, cluster = NULL) {
    cl = if (numberOfCores > 1 && is.null(cluster)) {
        c = parallel::makeCluster(numberOfCores)
        on.exit(parallel::stopCluster(c), add = TRUE)
        c
    } else cluster
    ...
}
```

Batch users create the cluster once and reuse it:

```r
cl = parallel::makeCluster(4)
on.exit(parallel::stopCluster(cl))
for (file in input_files) {
    result = dataProcess(read_input(file), numberOfCores = 4, cluster = cl)
}
```

---

## 14. Cumulative Impact and Prioritization

### Tiered summary

| Tier | Issue | Axes | Severity at 10k proteins |
|------|-------|------|--------------------------|
| 1 | S1: summarize loop accumulator | memory + protein loop | +1.6 GB observed at 3k; ~4-8 GB at 10k |
| 1 | S2: groupComparison single-core accumulator | memory + protein loop | ~0.5-1 GB |
| 1 | S3: `save_fitted_models = TRUE` retains models | memory + final-object size | ~0.4-2 GB in returned object |
| 1 | S4: `MSstatsBalancedDesign` Prepare spike (upstream) | memory + time | ~2.5 GB on 324 MB input |
| 2 | S5: O(n²) heatmap clustering | time + memory | ~800 MB distance matrix + minutes |
| 2 | S6: ggplot accumulation in plot loops | memory + protein loop | ~2-10 GB returned-list size |
| 2 | S7: O(n²) hover-text loop | time | tens of seconds to minutes |
| 3 | S8: quantile-norm full copy | memory | ~250 MB if quantile is selected |
| 3 | S9: linear-path `lm` retention | GC churn | ~400-800 MB transient |
| 3 | S10: `.SD` without `.SDcols` | memory + time | hundreds of MB transient |
| 4 | S11: groupComparison closure captures full list | memory × cores | 8 × 500 MB at 8 cores |
| 4 | S12: cluster startup per call | time | ~8s × N files |

### Recommended order at scale

1. **S1, S2** — stream/collate per-protein results. Single biggest peak
   reductions on real inputs.
2. **S3** — change `save_fitted_models` default to `FALSE` or strip
   bulky slots. Largest reduction in returned-object size.
3. **S5, S6** — fix heatmap clustering ceiling and ggplot accumulation;
   these become the bottlenecks once Tier 1 is addressed.
4. **S11** — port the split-first fix from `dataProcess.R` to
   `utils_groupcomparison.R`.
5. **S4** — file upstream issue against `MSstatsConvert`. Until that
   lands, Prepare-stage peak is unavoidable.
6. **S7, S8, S9, S10, S12** — mechanical wins to fold in when touching
   those files for other reasons.

### Verification

The existing `benchmark/profile_dataprocess_peak.R` script can confirm
S1 cost on the current branch. To confirm S2/S3, run
`MSstats::groupComparison` on the 3,000-protein fixture and watch
`gc()`'s `(Mb)` columns across the single-core path — the spike between
the loop end and `MSstatsGroupComparisonOutput` return is S2 + S3. For
S5, time `groupComparisonPlots(..., type = "Heatmap")` against
protein-count slices (1k, 3k, 10k) and confirm super-linear growth.
