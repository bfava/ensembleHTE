# AGENTS.md

`ensembleHTE` is an R package for ensemble-based heterogeneous treatment
effect (HTE) estimation and prediction with valid downstream inference,
using repeated K-fold cross-fitting. See
[README.md](https://bfava.com/ensembleHTE/README.md) for overview/quick
start and
[CONTRIBUTING.md](https://bfava.com/ensembleHTE/CONTRIBUTING.md) for
general contribution rules.

## Architecture

Two user-facing entry points drive everything:

- [R/ensemble_hte.R](https://bfava.com/ensembleHTE/R/ensemble_hte.R) —
  [`ensemble_hte()`](https://bfava.com/ensembleHTE/reference/ensemble_hte.md):
  causal/HTE estimation (R-/T-/S-/X-learner).
- [R/ensemble_pred.R](https://bfava.com/ensembleHTE/R/ensemble_pred.R) —
  [`ensemble_pred()`](https://bfava.com/ensembleHTE/reference/ensemble_pred.md):
  plain prediction with cross-fitting.

Both produce S3 objects (`ensemble_hte_fit`, `ensemble_pred_fit`)
consumed by analysis functions:

- [R/analysis_hte.R](https://bfava.com/ensembleHTE/R/analysis_hte.R) —
  [`blp()`](https://bfava.com/ensembleHTE/reference/blp.md),
  [`gates()`](https://bfava.com/ensembleHTE/reference/gates.md),
  [`clan()`](https://bfava.com/ensembleHTE/reference/clan.md)
- [R/analysis_pred.R](https://bfava.com/ensembleHTE/R/analysis_pred.R) —
  [`blp_pred()`](https://bfava.com/ensembleHTE/reference/blp_pred.md),
  [`gavs()`](https://bfava.com/ensembleHTE/reference/gavs.md)
- [R/compare_restricted.R](https://bfava.com/ensembleHTE/R/compare_restricted.R)
  —
  [`gates_restricted()`](https://bfava.com/ensembleHTE/reference/gates_restricted.md),
  [`gavs_restricted()`](https://bfava.com/ensembleHTE/reference/gavs_restricted.md)

Supporting modules:

- [R/ml_training.R](https://bfava.com/ensembleHTE/R/ml_training.R) —
  dispatcher to `mlr3` (most algorithms) or `grf` (causal forests);
  handles
  [`.encode_factors_if_needed()`](https://bfava.com/ensembleHTE/reference/dot-encode_factors_if_needed.md)
  via `fastDummies`.
- [R/utils.R](https://bfava.com/ensembleHTE/R/utils.R) —
  `create_folds()` (supports stratification + cluster/panel IDs),
  grouping, SE helpers.
- [R/news.R](https://bfava.com/ensembleHTE/R/news.R) —
  [`ensemble_news()`](https://bfava.com/ensembleHTE/reference/ensemble_news.md)
  /
  [`.fetch_update_status()`](https://bfava.com/ensembleHTE/reference/dot-fetch_update_status.md)
  reads
  [inst/update_status.dcf](https://bfava.com/ensembleHTE/inst/update_status.dcf).

S3 methods (`print.*`, `plot.*`, `summary.*`) are exported via roxygen —
see [NAMESPACE](https://bfava.com/ensembleHTE/NAMESPACE).

## Conventions

- **Style**: tidyverse style guide; roxygen2 for all docs
  (`@keywords internal` for non-exported helpers).
- **Internal helpers** are prefixed with `.` (e.g. `.fit_model_grf`,
  `.blp_single`) and live alongside the exported function that calls
  them.
- **Algorithm names** passed via `algorithms = c(...)`: `"grf"` routes
  to `grf`; everything else (`"lm"`, `"glmnet"`, `"xgboost"`,
  `"ranger"`, `"gbm"`, …) routes through `mlr3learners`.
- **Cross-fitting parameters** are always `M` (repetitions) × `K`
  (folds); per-repetition estimates are averaged with SEs combined
  across splits per Fava (2025).
- **Suggested packages** (`glmnet`, `gbm`, `ranger`) must be guarded
  with [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) or
  `skip_if_not_installed()` in tests.
- **Datasets**: `microcredit` is `LazyData`; rebuild via
  [data-raw/microcredit.R](https://bfava.com/ensembleHTE/data-raw/microcredit.R).

## Build / Test / Check

Run from the package root:

``` r

devtools::load_all(".")        # iterative development
devtools::document()           # regenerate man/*.Rd + NAMESPACE after roxygen edits
devtools::test()               # run testthat suite
devtools::check()              # full R CMD check — must be 0/0/0 (see cran-comments.md)
```

Heavy tests use `skip_on_cran()` and small `M`/`K` (e.g. `M=2, K=2`) —
keep new tests fast and CRAN-safe. Tests live in
[tests/testthat/](https://bfava.com/ensembleHTE/tests/testthat/) and are
organized by module (`test-ensemble_hte.R`, `test-analysis_functions.R`,
`test-validation.R`, …).

## Vignettes (expensive)

Vignettes do **not** train models at build time. They load pre-computed
`.rds` fits from
[vignettes/precomputed/](https://bfava.com/ensembleHTE/vignettes/precomputed/).
After changing the public API of
[`ensemble_hte()`](https://bfava.com/ensembleHTE/reference/ensemble_hte.md)
/
[`ensemble_pred()`](https://bfava.com/ensembleHTE/reference/ensemble_pred.md)
or the `microcredit` data, regenerate them:

``` r

source("data-raw/precompute_vignette_fits.R")
```

[data-raw/precompute_walkthrough.R](https://bfava.com/ensembleHTE/data-raw/precompute_walkthrough.R)
is the analog for the complete-guide article.

## Pitfalls

- Small `n` × small `K` with R-learner can yield zero-variance residuals
  — tests use `metalearner = "t"` to dodge this.
- `prop_score` must be supplied (or estimable) for
  [`ensemble_hte()`](https://bfava.com/ensembleHTE/reference/ensemble_hte.md);
  clipping is the caller’s responsibility.
- After editing any `@param`/`@return`/`@export`, run
  `devtools::document()` before committing — `man/*.Rd` and `NAMESPACE`
  are generated, never hand-edited.
- The local `ensembleHTE.Rcheck/` directory is `R CMD check` output;
  never commit edits to it.
- `.DS_Store` files appear under [R/](https://bfava.com/ensembleHTE/R/)
  and [vignettes/](https://bfava.com/ensembleHTE/vignettes/) — ignore,
  don’t track.

## Open work

Tracked in [TODO.md](https://bfava.com/ensembleHTE/TODO.md).
