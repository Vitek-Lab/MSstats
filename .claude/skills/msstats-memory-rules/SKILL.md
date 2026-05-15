---
name: msstats-memory-rules
description: Use when editing the MSstats R codebase — especially any change that touches a feature-level data.table, fits lm/survreg/lmerMod in a per-protein loop, or modifies the parallel summarize/groupComparison path. Surfaces the 12 anti-patterns extracted from real bugs on the Fix-Memory-Management branch — patterns that look correct, pass tests, and silently allocate hundreds of MB on real-scale inputs.
---

# MSstats memory & scalability rules

Twelve anti-patterns extracted from the bug fixes on the
`Fix-Memory-Management` branch. Each looked correct, passed tests, and
silently allocated hundreds of MB to multiple GB on real-scale inputs
(10k+ proteins, multi-GB feature tables).

**Full Wrong/Right code examples + commit hashes live in
`docs/CRITICAL-RULES.md` in the repo root.** Read that file for any rule
your edit touches.

## Rule index

| # | Rule (one-line) |
|---|-----------------|
| 1 | Drop columns by reference (`data.table::set(dt, j=, value=NULL)`), never by `dt[, !cols, with=FALSE]` |
| 2 | Update columns by keyed match (`other[input, on=, which=TRUE]` + `set()`), never by `merge()` on the full table |
| 3 | Use subset-write `dt[cond, col := val]`, never `dt[, col := ifelse(cond, val, col)]` |
| 4 | Convert to data.frame at output via `data.table::setDF(dt)`, never `as.data.frame(dt)` |
| 5 | Add columns via `df[["col"]] <- ...`, never `df <- data.frame(df, newcol = ...)` (especially in loops) |
| 6 | Sort with `data.table::setorder(dt, ...)`, never `dt <- dt[order(...), ]` |
| 7 | `rm()` fitted-model objects (lm/survreg) immediately after last use; strip `$model`, `$qr`, `$y` in the fitter |
| 8 | `rm()` function-parameter variables after their last read (`raw`, `peptides_dict`) |
| 9 | Combine column-subset + transform into a single `dt[, list(...)]` expression, never two passes |
| 10 | Match on the key tuple via `on=c("FEATURE", "FRACTION")`, never synthesise a `paste()`-concatenated key column |
| 11 | Null out consumed list slots (`list[[i]][[k]] <- NULL`) before downstream allocations |
| 12 | Pass per-element data via `parLapply(cl, X, FUN)`'s `X` argument; never reference large variables by closure |

## Pre-edit diagnostic checklist

Run through this before writing any patch that touches `R/`:

1. **Did I write `dt = dt[, ...]` anywhere?** Is it dropping columns
   (Rule 1), sorting (Rule 6), or column-subsetting (Rule 9)? Each has
   an in-place replacement.
2. **Did I call `merge()` on `input` or another feature-level table?**
   Can it be a keyed match instead (Rule 2)?
3. **Did I write `ifelse()` inside `[, := ...]`?** Almost always wrong
   (Rule 3) — `ifelse` allocates a full-length vector regardless of
   selectivity.
4. **Did I write `as.data.frame()` or `data.frame(dt, ...)`?** Use
   `setDF()` or `[[<-` instead (Rules 4, 5).
5. **Does a new function parameter remain in scope long after its last
   read?** `rm()` it (Rule 8). R cannot GC bound arguments.
6. **Did I introduce a fitted model (`lm`, `survreg`, `lmerMod`) in a
   per-protein loop?** Strip bulky slots in the fitter, `rm()` it right
   after `predict()` / `vcov()` (Rule 7).
7. **Did I add a `parLapply` or `clusterExport`?** Make sure the worker
   function doesn't reference a large variable by closure (Rule 12).
8. **Does my change synthesise a `paste()` column for matching?** Match
   on the tuple directly via `on=c(...)` instead (Rule 10).

## When in doubt — verify with profmem

If you're unsure whether a `[, := ...]` line is doing what you think:

```r
profmem::profmem({
    input[is.na(newABUNDANCE), nonmissing_orig := TRUE]
})
```

`profmem` reports every allocation > 0 bytes. On a multi-million-row
input, any single line that allocates more than a few MB is almost
certainly violating one of the rules above.

For stage-by-stage peak comparison across the whole pipeline:
`benchmark/profile_dataprocess_peak.R`.

## Why these rules exist

These aren't style preferences. Each rule maps to a specific bug fixed
on this branch — see `docs/CRITICAL-RULES.md` for the commit hash and
the memory cost that prompted the fix. The smallest of these saved
~40 MB per call; the largest cut peak memory by over 1 GB.

If you find a violation that isn't covered, add it to
`docs/CRITICAL-RULES.md` and update this skill's rule index.
