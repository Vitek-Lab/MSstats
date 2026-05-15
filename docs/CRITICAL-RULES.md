# CRITICAL RULES — Memory & Scalability Anti-Patterns

Rules extracted from real bugs fixed on the `Fix-Memory-Management`
branch. Each rule is a pattern that **looked correct, passed tests, and
silently allocated hundreds of MB** on real-scale inputs. The "Wrong"
column is what the original code did; the "Right" column is what
replaced it; the commit hash is where the lesson was learned.

If you are editing this codebase and your change touches a feature-level
`data.table`, scan this file before writing the patch. Most of these
rules are violated easily and are not caught by tests.

---

## Rule index

| # | Rule | Cost when violated |
|---|------|--------------------|
| 1 | Never drop columns with `dt[, !cols, with = FALSE]` | Full deep-copy of every retained column |
| 2 | Never use `merge()` on the full feature table to update a column | 2-3 simultaneous copies of `input` |
| 3 | Never use `ifelse()` for conditional `:=` updates | Allocates a full-length vector regardless of how many rows match |
| 4 | Never use `as.data.frame()` on a large `data.table` at output | Deep-copies every column |
| 5 | Never use `data.frame(input, newcol = ...)` to add columns inside a loop | One full copy per call |
| 6 | Never sort by `dt[order(...), ]` | Builds a full row-reshuffled copy |
| 7 | Never let `lm` / `survreg` outlive their last use in a per-protein loop | Model object retains `$model`, `$qr`, `$residuals` |
| 8 | Never leave function-parameter variables alive after their last read | `raw`, `peptides_dict` sit ~300 MB through entire summarization |
| 9 | Never do column-subset + transform as two steps | Two new data.tables; combine into one expression |
| 10 | Never synthesise a `paste()` key column just to filter/match | Two large character vectors + paste cost; match on the tuple directly |
| 11 | Never hold a nested per-protein list after extracting from it | Per-protein survival data.tables sit alive through `.finalizeInput` |
| 12 | Never capture a large variable via `parLapply` closure | Every worker receives a full copy by serialization |

---

## Rule 1. Drop columns by reference, not by subsetting

**Wrong.** Builds a brand-new `data.table` and **deep-copies every
retained column**.

```r
input = input[, !(colnames(input) %in% c("ABUNDANCE_RUN", "ABUNDANCE_FRACTION")),
              with = FALSE]
```

**Right.** Removes the column pointer in place; zero allocation.

```r
data.table::set(input, j = "ABUNDANCE_RUN", value = NULL)
data.table::set(input, j = "ABUNDANCE_FRACTION", value = NULL)
```

For many columns:

```r
drop_cols = setdiff(colnames(input), output_cols)
for (col in drop_cols) data.table::set(input, j = col, value = NULL)
```

**Where this bit us.** `b9e89d4`, `12a0380`, `c6cc9ae`. Four normalize
sites + the prepare stage + the output stage were each doing a ~250 MB
transient copy to drop 1-2 temporary columns.

---

## Rule 2. Update columns by keyed match, not by `merge()`

**Wrong.** `merge()` allocates a new `data.table`. The column-subset
inside compounds it (Rule 1).

```r
input = merge(input[, colnames(input) != "newABUNDANCE", with = FALSE],
              predicted_survival,
              by = setdiff(cols, "newABUNDANCE"),
              all.x = TRUE)
```

**Right.** Look up the match indices once; write the values in place.

```r
idx = predicted_survival[input, on = join_cols, which = TRUE,
                         mult = "first"]
data.table::set(input, j = "newABUNDANCE",
                value = predicted_survival$newABUNDANCE[idx])
data.table::set(input, j = "predicted",
                value = predicted_survival$predicted[idx])
```

**Where this bit us.** `745f62b`, `0e83ce1`, `1d2dfd9`. The
`.finalizeTMP` merge alone was triggering ~1.1 GB of simultaneous copies
on a 350 MB input. Note that even `input[other, col := i.col, on = ...]`
("update by join") allocates table-sized vectors internally — prefer
match-indices + `set()` for the largest tables.

---

## Rule 3. Subset-write with `:=`, never `ifelse()`

**Wrong.** `ifelse()` allocates a full-length result vector for every
row, regardless of how many actually match the condition.

```r
input[, nonmissing_orig := ifelse(is.na(newABUNDANCE), TRUE, nonmissing_orig)]

single_protein[, predicted := ifelse(censored & (LABEL == "L"),
                                     predicted, NA)]
```

**Right.** Subset-targeted `:=` only touches rows that change.

```r
input[is.na(newABUNDANCE), nonmissing_orig := TRUE]

single_protein[!(censored & LABEL == "L"), predicted := NA]
single_protein[censored & LABEL == "L",
               newABUNDANCE := predicted]
```

**Where this bit us.** `f55485e`, `53a8f14`. Runs once per protein × per
finalize pass — looks like nothing in isolation, but the allocations
compound across the summarize loop into measurable peak.

---

## Rule 4. Use `setDF()`, not `as.data.frame()`, when returning

**Wrong.** Deep-copies every column vector.

```r
list(FeatureLevelData = as.data.frame(input),
     ProteinLevelData = as.data.frame(rqall),
     SummaryMethod = method)
```

**Right.** Changes the class in place. Zero allocation.

```r
data.table::setDF(input)
data.table::setDF(rqall)
list(FeatureLevelData = input,
     ProteinLevelData = rqall,
     SummaryMethod = method)
```

`data.table` already inherits from `data.frame`, so consumers that test
`is.data.frame()` will be fine either way.

**Where this bit us.** `3197f4e`. The output step was momentarily
holding two full copies (the data.table and its as.data.frame twin) for
data worth ~400 MB.

---

## Rule 5. Add columns by reference, not by `data.frame(...)`

**Wrong.** `data.frame(input, newcol = ...)` builds a new data.frame by
copying every existing column.

```r
input = data.frame(input, "abs.resids" = abs.resids, "fitted" = fitted)
input = data.frame(input, "loess.fitted" = loess.fitted)
```

**Right.** Convert once, then assign columns in place.

```r
input = as.data.frame(input)                       # one conversion
input[["abs.resids"]]   = abs(fit$residuals)
input[["fitted"]]       = fit$fitted.values
input[["loess.fitted"]] = fitted(fit.loess)
input[["abs.resids"]]   = NULL                     # drop in place too
```

(If you need data.table semantics throughout, stay in data.table and use
`:=` instead — but **do not mix**. `data.frame(dt, ...)` silently
strips data.table-ness.)

**Where this bit us.** `6a51032`. The unequal-variance loop was doing
~5 full-copy operations per iteration per protein.

---

## Rule 6. Sort with `setorder()`, not `dt[order(...), ]`

**Wrong.** Builds a full row-reshuffled copy of the table.

```r
input = input[order(LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL,
                    RUN, PROTEIN, PEPTIDE, TRANSITION), ]
```

**Right.** Rewrites the column vectors in place.

```r
data.table::setorder(input, LABEL, GROUP_ORIGINAL, SUBJECT_ORIGINAL,
                     RUN, PROTEIN, PEPTIDE, TRANSITION)
```

**Where this bit us.** `c6cc9ae`. Prepare stage was allocating an entire
input-sized copy just to reorder rows.

---

## Rule 7. Free fitted-model objects immediately after their last use

`lm` and `survreg` return objects that retain the **full fitting data**:

- `lm$model` — the data frame used for fitting (all formula columns)
- `lm$qr$qr`, `lm$residuals`, `lm$fitted.values`, `lm$effects`
- `survreg$y`, `survreg$linear.predictors`
- The formula environment — can chain to arbitrary enclosing scope

In a per-protein loop, every kilobyte you hold past the last `predict()`
call multiplies by the number of proteins and increases GC churn.

**Wrong.**

```r
survival_fit = .fitSurvival(...)
predict(survival_fit, newdata = .SD)
# ... 40 more lines of protein body ...
return(list(result, survival))    # survival_fit only released here
```

**Right.** Two complementary moves:

```r
survival_fit = .fitSurvival(...)
predict(survival_fit, newdata = .SD)
rm(survival_fit)                  # release immediately
```

Inside `.fitSurvival` / `.fitLinearModel`, strip the bulky slots from
what you return:

```r
fit$y                = NULL
fit$linear.predictors = NULL
# or for lm:
fit$model            = NULL
```

Only `coef(fit)`, `vcov(fit)`, `fit$rank`, `fit$df.residual` are
typically needed downstream.

**Where this bit us.** `420a380`. The TMP path was patched; the linear
path still holds `lm` across the protein body — see S9 in
`scalability_issues.md`.

---

## Rule 8. `rm()` function-parameter variables after their last read

R cannot garbage-collect a function argument as long as it is bound in
the active call frame. If `dataProcess(raw, ...)` reads `raw` once at
line 146 but returns at line 177, `raw` (~250-400 MB) sits alive through
the entire summarization step.

**Wrong.**

```r
dataProcess = function(raw, ...) {
    input = MSstatsPrepareForDataProcess(raw, logTrans, fix_missing)
    # raw never read again, but stays alive ~350 MB
    input = MSstatsNormalize(input, normalization, peptides_dict, nameStandards)
    # peptides_dict never read again
    ...
    # summarization runs with raw + peptides_dict both still resident
}
```

**Right.**

```r
dataProcess = function(raw, ...) {
    input = MSstatsPrepareForDataProcess(raw, logTrans, fix_missing)
    rm(raw)
    input = MSstatsNormalize(input, normalization, peptides_dict, nameStandards)
    rm(peptides_dict)
    ...
}
```

**Where this bit us.** `98950b4`.

---

## Rule 9. Combine column-subset + transform into one expression

**Wrong.** Two passes, two new data.tables.

```r
cols  = intersect(c("PROTEIN", "PEPTIDE", "FEATURE", ...), colnames(input))
input = input[, cols, with = FALSE]                      # COPY 1: ~120 MB
data.table::setnames(input, "censored", "is_censored")
input = input[, list(protein = as.character(PROTEIN),
                     peptide = as.character(PEPTIDE),
                     ...)]                               # COPY 2: ~100 MB
```

**Right.** Build the working copy in a single `data.table` expression.

```r
has_censored = is.element("censored", colnames(input))
work = input[, list(protein  = as.character(PROTEIN),
                    peptide  = as.character(PEPTIDE),
                    feature  = as.character(FEATURE),
                    run      = as.character(originalRUN),
                    label    = as.character(LABEL),
                    log2inty = ifelse(!(is.na(ABUNDANCE) |
                                       if (has_censored) censored else FALSE),
                                      ABUNDANCE, NA),
                    is_obs   = FALSE)]
work[, is_obs := !is.na(log2inty)]
```

**Where this bit us.** `0a1b515` (`.selectHighQualityFeatures`).

---

## Rule 10. Match on the key tuple, never on a `paste()`-concatenated key

**Wrong.** Synthesizes two large character vectors and depends on a
substring of one matching a substring of the other.

```r
select_fraction[, tmp := paste(FEATURE, FRACTION, sep = "_")]
input$tmp = paste(input$FEATURE, input$FRACTION, sep = "_")
input = input[tmp %in% select_fraction$tmp, ]
```

**Right.** Match on the pair of columns directly. No string allocation.

```r
keep_idx = select_fraction[input,
                           on = c("FEATURE", "FRACTION"),
                           which = TRUE, mult = "first"]
input = input[!is.na(keep_idx)]
```

**Where this bit us.** `1d2dfd9`. The two `paste()` calls alone were
allocating two ~200k-string character vectors per call, and the
character-key compare was much slower than the integer-key match.

---

## Rule 11. Null out consumed list slots before downstream work

When you have a nested list `list(list(summary_i, survival_i), ...)` and
the next stage only needs the survivals collated, drop the references
to the per-protein survivals **before** doing more allocations.

**Wrong.** All three of `summarized`, `protein_summaries`, and
`predicted_survival` coexist across `.finalizeInput`.

```r
predicted_survival = data.table::rbindlist(lapply(summarized, function(x) x[[2]]))
protein_summaries  = lapply(summarized, function(x) x[[1]])
input = .finalizeInput(input, predicted_survival, ...)
```

**Right.** Null out the consumed slot before the next allocation.

```r
predicted_survival = data.table::rbindlist(lapply(summarized, function(x) x[[2]]))
for (i in seq_along(summarized)) summarized[[i]][[2]] = NULL    # release survivals
protein_summaries = lapply(summarized, function(x) x[[1]])
rm(summarized)
input = .finalizeInput(input, predicted_survival, ...)
```

**Where this bit us.** `d622fe4`, `9302f25`.

---

## Rule 12. Pass per-element data through `parLapply(X, ...)`, never via closure

`parallel::parLapply(cl, X, FUN)` serializes `FUN` **including its
enclosing environment**. Any large variable that `FUN` references by
name — even if you didn't `clusterExport` it — gets sent to every worker.

**Wrong.** `summarized_list` is captured by closure and copied to every
worker.

```r
test_results = parallel::parLapply(cl, all_proteins_id, function(i) {
    MSstatsGroupComparisonSingleProtein(i, summarized_list[[i]], ...)
})
```

**Right.** Pre-split, then iterate over the slices.

```r
protein_data_list = lapply(all_proteins_id, function(i) summarized_list[[i]])
test_results = parallel::parLapply(cl, protein_data_list,
    function(single_protein_data) {
        MSstatsGroupComparisonSingleProtein(single_protein_data, ...)
    })
```

**Where this bit us.** The `dataProcess` parallel path was fixed this
way (commit `98950b4`-era refactor); the same pattern still exists in
`utils_groupcomparison.R:475-491` — see S11 in `scalability_issues.md`.

---

## Diagnostic checklist for new edits

Before merging any change that touches a feature-level `data.table`:

1. **Did I write `dt = dt[, ...]` anywhere?** If so — is it dropping
   columns (Rule 1), sorting (Rule 6), or subsetting (Rule 9)? Each has
   an in-place replacement.
2. **Did I call `merge()` on `input`?** Can it be a keyed match instead
   (Rule 2)?
3. **Did I write `ifelse()` inside `[, := ...]`?** Almost always wrong
   (Rule 3).
4. **Did I write `as.data.frame()` or `data.frame(dt, ...)`?** Use
   `setDF()` or `[[<-` (Rules 4, 5).
5. **Does a new function parameter remain in scope long after its last
   read?** `rm()` it (Rule 8).
6. **Did I introduce a fitted model (`lm`, `survreg`, `lmerMod`) in a
   per-protein loop?** Strip bulky slots in the fitter, `rm()` it right
   after `predict()` / `vcov()` (Rule 7).
7. **Did I add a `parLapply` or `clusterExport`?** Make sure the worker
   function doesn't reference a large variable by closure (Rule 12).
8. **Does my change synthesise a `paste()` column?** Match on the tuple
   instead (Rule 10).

## Verifying suspected violations

Use `profmem::profmem()` around the suspect block — it reports every
allocation `> 0 bytes`. For peak memory, the
`benchmark/profile_dataprocess_peak.R` script in this repo gives
stage-by-stage `gc()` deltas across the whole pipeline.

If a single `:=` or `[, ...]` line allocates more than a few MB on a
multi-million-row input, it is violating one of the rules above.
