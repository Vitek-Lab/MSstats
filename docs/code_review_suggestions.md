# Code review suggestions — MSstats branch `Fix-Memory-Management`

Findings from the debugging + memory-tuning session on this branch. Each
item follows the same shape: **the issue**, **the effect**, **a solution**.

Sections: code quality, memory, time.

---

## Code quality

### 1. Join key list is duplicated in three places

**The issue.** The same list of "columns that uniquely identify a row in
per-protein survival" is written out three separate times:

- `R/dataProcess.R::MSstatsSummarizeSingleTMP` line 483
- `R/dataProcess.R::MSstatsSummarizeSingleLinear` line 377
- `R/utils_output.R::.finalizeTMP` line 146

Each list is supposed to agree with the other two. Nothing enforces that.

**The effect.** During this session, the three lists got out of sync
three different times — `ref` vs `ref_covariate`, missing `LABEL`,
missing `PROTEIN`. Each time, `dataProcess` crashed: once with a
"missing columns" merge error, once with a Cartesian-explosion error.
Same root cause, different symptoms, hours lost diagnosing.

**A solution.** Move the canonical list into one constant and have all
three sites read it.

```r
# R/utils_output.R (or a new R/utils_schema.R)
.SURVIVAL_JOIN_KEYS <- c("PROTEIN", "FEATURE", "RUN", "LABEL",
                         "ref_covariate", "cen")

.intersect_survival_keys <- function(available_cols) {
    intersect(available_cols, .SURVIVAL_JOIN_KEYS)
}
```

In all three sites:
```r
cols <- .intersect_survival_keys(colnames(single_protein))
```

One source of truth. Adding a new key is a one-line edit, not three.

---

### 2. Dead conditional logic in `.finalizeInput`

**The issue.** `R/utils_output.R::.finalizeInput` has a commented-out
`if (method == "TMP") … else .finalizeLinear(...)` block. The function
always calls `.finalizeTMP`, regardless of method. `.finalizeLinear`
exists but is never used.

**The effect.** Readers can't tell whether linear was deliberately
unified with TMP or whether the dispatch was accidentally broken.
`.finalizeLinear` is maintenance overhead with no caller.

**A solution.** Pick one path. Either:
- Delete `.finalizeLinear` and the commented block, and add a one-line
  comment explaining that both methods share the TMP finalization path.
- Or restore the dispatch and audit what's different between the two
  finalize functions.

Dead code is worse than removed code — it suggests behaviour that
doesn't exist.

---

### 3. Pipeline stages are duplicated in `dataProcess` and the profile script

**The issue.** The benchmark script
`benchmark/profile_dataprocess_peak.R` replays the same stage sequence
as `dataProcess` — `MSstatsPrepareForDataProcess → MSstatsNormalize →
MSstatsMergeFractions → …`. Two copies, kept in sync by hand.

**The effect.** If `dataProcess` reorders or adds a stage, the profile
script silently profiles the wrong pipeline and reports misleading
numbers.

**A solution.** Add a callback hook to `dataProcess`:

```r
dataProcess <- function(raw, ..., .stage_callback = NULL) {
    record <- function(label) if (!is.null(.stage_callback)) .stage_callback(label)
    ...
    input <- MSstatsPrepareForDataProcess(raw, logTrans, fix_missing)
    record("after Prepare")
    ...
}
```

Then the profile script becomes a few lines using the real `dataProcess`
pipeline — no duplication.

---

### 4. Linter output is mostly noise

**The issue.** Every file edit emits dozens of style warnings — `=` vs
`<-`, 4-space vs 2-space indentation, line length, snake_case. The
package uses 4-space indent and `=` consistently, so the warnings are
fighting the established style.

**The effect.** Real warnings (no-visible-binding, unused-variable,
missing imports) get buried in the cosmetic noise. People stop reading
lintr output entirely.

**A solution.** Decide on a house style once, then commit a `.lintr`
file that disables the rules you don't care about:

```yaml
linters: with_defaults(
    assignment_linter      = NULL,
    indentation_linter     = NULL,
    line_length_linter     = NULL,
    object_name_linter     = NULL,
    object_usage_linter    = lintr::object_usage_linter
)
```

The warnings that remain are the ones worth fixing.

---

### 5. Parallel workers run a different version of MSstats than the parent

**The issue.** Running tests under `pkgload::load_all(".")` loads dev
code into the parent R session. `dataProcess(..., numberOfCores = 2)`
spawns worker processes that don't see the dev source — they load
whatever's installed in the R library.

**The effect.** Workers run an older version of the per-protein
summarize functions. Their outputs disagree with what the parent's
`.finalizeTMP` expects. Three crashes in this session traced back to
this. The fix every time was "remember to `R CMD INSTALL .`", which is
easy to forget.

**A solution.** In `MSstatsSummarizeWithMultipleCores`, after
`makeCluster`, detect dev mode and sync the workers:

```r
cl <- parallel::makeCluster(numberOfCores)
if (requireNamespace("pkgload", quietly = TRUE) &&
    pkgload::is_dev_package("MSstats")) {
    pkg_path <- pkgload::pkg_path()
    parallel::clusterCall(cl, function(p) {
        pkgload::load_all(p, quiet = TRUE)
    }, pkg_path)
}
```

About five lines. Workers automatically use the same code as the parent
during development. No effect in production (the `is_dev_package`
check is false when MSstats is installed normally).

---

### 6. Tests assert on `ncol(...)` instead of on required columns

**The issue.** Several tests check the exact number of columns in
intermediate tables — e.g. `expect_equal(ncol(protein_level_summary),
4)`.

**The effect.** Every time the schema legitimately adds a column (like
`LABEL` or `ref_covariate`), these tests fail with a number mismatch
even though no real contract was violated. Test maintenance after a
schema change becomes a chore.

**A solution.** Assert on the contract instead — which columns must be
present.

```r
required_cols <- c("Protein", "RUN", "LogIntensities", "Variance")
expect_true(all(required_cols %in% colnames(protein_level_summary)),
            info = "ProteinLevelData must contain required columns")
```

Adding a column doesn't break the test. Removing a required column
still does.

For the snapshot equality test, restrict the comparison to a stable
column subset:

```r
stable_cols <- c("PROTEIN", "FEATURE", "RUN", "ABUNDANCE", "newABUNDANCE")
expect_true(data.table::fsetequal(dt1[, ..stable_cols], dt2[, ..stable_cols]))
```

---

## Memory

### 7. `.SD` without `.SDcols` materialises every column

**The issue.** Several places use `.SD` (the "subset of data" handle in
`data.table[, := , by = ...]`) without restricting which columns it
should expose. When the callee only reads a handful of columns, `.SD`
hands it all 25 — for no reason.

Affected sites:
- `R/utils_feature_selection.R:148` and `:202`
- `R/utils_normalize.R:54` and `:153`
- `R/utils_output.R:163` and `:198`
- `R/utils_summarization_prepare.R:128`

**The effect.** Each call materialises a list of every column for the
grouping. On a 1 MB-per-column table with 20 columns and a callee that
only reads 6, that's 14 columns × number of groups of wasted allocation
per call. Small per call, adds up across the pipeline on large inputs.

**A solution.** Look at what each callee actually reads. Pass exactly
those columns via `.SDcols`.

```r
input[, NonMissingStats := .getNonMissingFilterStats(.SD, censored_symbol),
      .SDcols = c("LABEL", "newABUNDANCE", "censored", "INTENSITY",
                  "n_obs_run", "n_obs")]
```

Now `.SD` is six columns wide, not twenty-five.

---

### 8. The summarize loop holds all per-protein results in memory before collation

**The issue.** `MSstatsSummarizeWithSingleCore` builds a list of length
`num_proteins` and assigns each iteration's result into one slot:

```r
summarized_results <- vector("list", num_proteins)
for (protein_id in seq_len(num_proteins)) {
    summarized_results[[protein_id]] <- MSstatsSummarizeSingleTMP(...)
}
```

Every iteration's `list(summary, survival)` stays alive until the loop
ends and the consumer runs `rbindlist`.

**The effect.** Memory grows linearly through the loop. On a 500-protein
fixture, `used` Vcells climbed from 219 to 417 MB across the loop. On a
3000-protein fixture, this drove a +1.6 GB spike — by far the biggest
remaining peak contributor at real-data scale.

**A solution.** Two options, depending on how much API change is
acceptable:

- **Big change (breaks API):** make `MSstatsSummarizeWithSingleCore`
  collate as it goes — accumulate two flat lists (summaries, survivals)
  and return them already `rbindlist`-ed. Bump the package's major
  version because the documented return shape changes.

- **Small change (no API break):** keep the nested-list return shape,
  but in `MSstatsSummarizationOutput`, null out
  `summarized[[i]][[2]]` immediately after the survival rbindlist
  consumes it. Halves the simultaneously-live survival data. (Some of
  this is already on the branch — verify it survived recent edits.)

The big change is the one that meaningfully lowers peak on real-scale
inputs.

---

### 9. Remaining column-subset copies on the full feature-level table

**The issue.** `dt[, cols, with = FALSE]` makes a new data.table by
deep-copying each retained column. A few of these still run on the
full-size feature table:

- `R/utils_normalize.R:118` — quantile-normalisation branch's
  `merge(input[, ..., with = FALSE], ...)`. Not exercised by the
  default test fixture, since the default is median normalisation.
- `R/utils_output.R:75-79` — two-step `unique(input[..., cols, with =
  FALSE])[, cols != "GROUP", with = FALSE]`. Two intermediates where
  one would do.
- `R/utils_summarization.R:97` — `as.matrix(wide[, features, with =
  FALSE])`. Per-protein, small.

**The effect.** Each copy is `nrow × ncol × 8 bytes`. On the test
fixture none of these touch peak (they sit downstream of larger
spikes). On real-scale inputs they're real allocations.

**A solution.** For `utils_normalize.R:118` (quantile branch): replace
the `merge` with the keyed-match pattern that's already in
`.finalizeTMP`. Add a quantile-normalisation test fixture first so any
change is verifiable.

For `utils_output.R:75-79`: combine the two subset operations into one.
The `unique()` step can take `cols` excluding `GROUP` directly.

For `utils_summarization.R:97`: the matrix is needed for
`median_polish_summary`, so the copy is unavoidable. Leave it.

---

### 10. `MSstatsConvert::MSstatsBalancedDesign` dominates the Prepare-stage spike

**The issue.** `R/utils_checks.R::.checkUnProcessedDataValidity` calls
`MSstatsConvert::MSstatsBalancedDesign` around line 175. That single
function caused a 108 MB spike on a 10.8 MB input and a 2.5 GB spike on
a 324 MB input — by far the largest single source of peak memory in
`dataProcess`.

**The effect.** Even with everything else optimised, this stage sets the
floor for total peak on any input. It pads the feature/run/label grid to
a balanced design, which means materialising the full cross-product
before filling NAs.

**A solution.** This function isn't in MSstats. It's in the
`MSstatsConvert` package. The structural fix has to live there:

- Stream the padding row-by-row instead of materialising the full grid.
- Or skip the padding when the input is already balanced (cheap check
  up front).

File an issue or PR against `MSstatsConvert` linking to the profile
output. Until that lands, the Prepare spike is unavoidable from inside
`MSstats`.

---

## Time

### 11. `paste()` on full feature-level columns is slow

**The issue.** Several places build identifier columns by concatenating
two existing columns with `paste()`:

- `R/utils_checks.R:208-211` — `PEPTIDE` and `TRANSITION`
- `R/utils_checks.R:289` — `FEATURE`

On a 6-million-row table, each `paste()` allocates 6 million character
strings, each going through R's global string-interning hash.

**The effect.** Pure CPU cost. Doesn't affect peak memory much but does
make `MSstatsPrepareForDataProcess` slower than it has to be.

**A solution.** Two cheap wins:

- `paste0()` is faster than `paste(..., sep = "_")` because it skips the
  printf-style formatting layer. Whenever the separator is a fixed
  literal, prefer `paste0`.
- `stringi::stri_c()` is 2-4× faster than `paste0` for large character
  vectors because it skips R's recycling/coercion logic:

```r
input$PEPTIDE <- stringi::stri_c(input$PEPTIDESEQUENCE, "_",
                                 input$PRECURSORCHARGE)
```

`stringi` is often already an indirect dependency (data.table pulls it
in). Confirm with `tools::package_dependencies("MSstats", recursive =
TRUE)`.

---

### 12. Factor calls scan their columns multiple times in `.makeFactorColumns`

**The issue.** `R/utils_checks.R::.makeFactorColumns` does eight
`factor(...)` calls in sequence, then a `setorder`. The last call on
`RUN` computes `unique(RUN)` twice:

```r
input[, RUN := factor(RUN, levels = unique(RUN),
                      labels = seq_along(unique(RUN)))]
```

`unique(RUN)` is a full pass over the column — and the call does it
twice for the same answer.

**The effect.** Each `unique()` is O(nrow). On 6 million rows that's a
few hundred milliseconds wasted per call.

**A solution.** Compute once, reuse:

```r
run_levels <- unique(RUN)
input[, RUN := factor(RUN, levels = run_levels,
                      labels = seq_along(run_levels))]
```

The other seven `factor()` calls don't have this exact issue but are
worth scanning — anywhere `unique(x)` shows up twice in the same
expression, hoist it.

---

### 13. Per-iteration cleanup inside the summarize loop is incomplete

**The issue.** Each iteration of `MSstatsSummarizeSingleTMP` (and the
linear counterpart) creates a `survival_fit` object that's typically the
last allocation alive in the iteration. The model object can be 100+ KB
and isn't released until function return.

Two related spots:
- `R/dataProcess.R::MSstatsSummarizeSingleTMP` line ~561 — `predict(survival_fit, newdata = .SD)`. `.SD` materialises the whole `single_protein` even though `predict` only reads predictor columns.
- The pre-set-NA-then-overwrite pattern around line 525.

**The effect.** A few MB per protein × 3000 proteins = ~30-60 MB of
extra peak on the summarize stage that doesn't need to exist.

**A solution.** Three small edits:

- Add `.SDcols = c("FEATURE", "RUN", "is_labeled_ref")` (or whichever
  predictor columns the survival model uses) to the `predict(...,
  newdata = .SD)` call.
- Insert `rm(survival_fit)` immediately after the last `predict()` call
  in each branch, before the next allocation.
- Replace the "set to NA, then conditionally overwrite" pattern with a
  single write of the right values.

---

### 14. Cluster startup happens on every `dataProcess` call

**The issue.** Each `dataProcess(..., numberOfCores = 2)` does:

```r
cl <- parallel::makeCluster(numberOfCores)
parallel::clusterExport(cl, c(...))
... parLapply ...
parallel::stopCluster(cl)
```

`makeCluster` spawns fresh R processes — 0.5 to 2 seconds of startup
per worker. For batch scripts that call `dataProcess` on many input
files, that cost is paid every file.

**The effect.** A 100-file batch wastes 50-200 seconds just on cluster
setup and teardown.

**A solution.** Accept an optional pre-made cluster:

```r
dataProcess <- function(raw, ..., numberOfCores = 1, cluster = NULL) {
    cl <- if (numberOfCores > 1 && is.null(cluster)) {
        c <- parallel::makeCluster(numberOfCores)
        on.exit(parallel::stopCluster(c), add = TRUE)
        c
    } else cluster
    ...
}
```

Batch users can then create the cluster once:

```r
cl <- parallel::makeCluster(4)
on.exit(parallel::stopCluster(cl))
for (file in input_files) {
    result <- dataProcess(read_input(file), numberOfCores = 4, cluster = cl)
    ...
}
```

---

### 15. Memory-test thresholds are too loose to catch real regressions

**The issue.** The current memory-test thresholds in
`inst/tinytest/test_dataProcess.R` are:

```
RETAINED_LIMIT    = 5
PEAK_LIMIT_TMP    = 75
PEAK_LIMIT_LINEAR = 75
```

The actual numbers on the 10.8 MB fixture are ~1.5x retained, ~12-15x
peak. The thresholds are set at 3-5× safety margin.

**The effect.** A future change could double peak memory (12x → 24x) and
the tests would silently pass — 24x is still under 75x. The tests
aren't doing their regression-detecting job.

**A solution.** Tighten the thresholds to closer to the actual baselines
with maybe 1.5-2× margin:

```r
RETAINED_LIMIT    <- 3    # current 1.5x; this allows 2x slack
PEAK_LIMIT_TMP    <- 25   # current ~12x; this allows 2x slack
PEAK_LIMIT_LINEAR <- 25   # current ~15x; this allows ~1.6x slack
```

Run the test five or ten times locally to confirm normal noise doesn't
flake the new thresholds. If CI runs on a different machine, calibrate
there separately.

---

## Recommended order

If you're triaging, here's a sensible rough order, easy and impactful
first:

1. **Item 1** (single source of truth for join keys) — small change,
   prevents bug class that hit three times.
2. **Item 5** (workers load dev source) — small change, prevents
   recurring "I edited the source but parallel crashes" surprise.
3. **Item 15** (tighten test thresholds) — small change, makes tests
   useful again.
4. **Item 7** (`.SDcols` everywhere) — mechanical, small win
   throughout.
5. **Item 14** (reusable cluster argument) — medium effort, big win
   for batch users.
6. **Item 8** (streaming summarize accumulator) — bigger refactor,
   biggest win on real-scale inputs. Save for a major version bump.
7. **Item 10** (file `MSstatsConvert` issue) — out-of-scope but worth
   documenting.

Everything else (2, 3, 4, 6, 9, 11, 12, 13) is small change to fold in
when touching those files for other reasons.
