# AGENTS.md

`ensembleHTE` is an R package for ensemble-based heterogeneous treatment effect (HTE) estimation and prediction with valid downstream inference, using repeated K-fold cross-fitting. See [README.md](README.md) for overview/quick start and [CONTRIBUTING.md](CONTRIBUTING.md) for general contribution rules.

## Architecture

Two user-facing entry points drive everything:

- [R/ensemble_hte.R](R/ensemble_hte.R) — `ensemble_hte()`: causal/HTE estimation (R-/T-/S-/X-learner).
- [R/ensemble_pred.R](R/ensemble_pred.R) — `ensemble_pred()`: plain prediction with cross-fitting.

Both produce S3 objects (`ensemble_hte_fit`, `ensemble_pred_fit`) consumed by analysis functions:

- [R/analysis_hte.R](R/analysis_hte.R) — `blp()`, `gates()`, `clan()`
- [R/analysis_pred.R](R/analysis_pred.R) — `blp_pred()`, `gavs()`
- [R/compare_restricted.R](R/compare_restricted.R) — `gates_restricted()`, `gavs_restricted()`

Supporting modules:

- [R/ml_training.R](R/ml_training.R) — dispatcher to `mlr3` (most algorithms) or `grf` (causal forests); handles `.encode_factors_if_needed()` via `fastDummies`.
- [R/utils.R](R/utils.R) — `create_folds()` (supports stratification + cluster/panel IDs), grouping, SE helpers.
- [R/news.R](R/news.R) — `ensemble_news()` / `.fetch_update_status()` reads [inst/update_status.dcf](inst/update_status.dcf).

S3 methods (`print.*`, `plot.*`, `summary.*`) are exported via roxygen — see [NAMESPACE](NAMESPACE).

## Conventions

- **Style**: tidyverse style guide; roxygen2 for all docs (`@keywords internal` for non-exported helpers).
- **Internal helpers** are prefixed with `.` (e.g. `.fit_model_grf`, `.blp_single`) and live alongside the exported function that calls them.
- **Algorithm names** passed via `algorithms = c(...)`: `"grf"` routes to `grf`; everything else (`"lm"`, `"glmnet"`, `"xgboost"`, `"ranger"`, `"gbm"`, ...) routes through `mlr3learners`.
- **Cross-fitting parameters** are always `M` (repetitions) × `K` (folds); per-repetition estimates are averaged with SEs combined across splits per Fava (2025).
- **Suggested packages** (`glmnet`, `gbm`, `ranger`) must be guarded with `requireNamespace()` or `skip_if_not_installed()` in tests.
- **Datasets**: `microcredit` is `LazyData`; rebuild via [data-raw/microcredit.R](data-raw/microcredit.R).

## Build / Test / Check

Run from the package root:

```r
devtools::load_all(".")        # iterative development
devtools::document()           # regenerate man/*.Rd + NAMESPACE after roxygen edits
devtools::test()               # run testthat suite
devtools::check()              # full R CMD check — must be 0/0/0 (see cran-comments.md)
```

Heavy tests use `skip_on_cran()` and small `M`/`K` (e.g. `M=2, K=2`) — keep new tests fast and CRAN-safe. Tests live in [tests/testthat/](tests/testthat/) and are organized by module (`test-ensemble_hte.R`, `test-analysis_functions.R`, `test-validation.R`, ...).

## Vignettes (expensive)

Vignettes do **not** train models at build time. They load pre-computed `.rds` fits from [vignettes/precomputed/](vignettes/precomputed/). After changing the public API of `ensemble_hte()` / `ensemble_pred()` or the `microcredit` data, regenerate them:

```r
source("data-raw/precompute_vignette_fits.R")
```

[data-raw/precompute_walkthrough.R](data-raw/precompute_walkthrough.R) is the analog for the complete-guide article.

## Pitfalls

- Small `n` × small `K` with R-learner can yield zero-variance residuals — tests use `metalearner = "t"` to dodge this.
- `prop_score` must be supplied (or estimable) for `ensemble_hte()`; clipping is the caller's responsibility.
- After editing any `@param`/`@return`/`@export`, run `devtools::document()` before committing — `man/*.Rd` and `NAMESPACE` are generated, never hand-edited.
- The local `ensembleHTE.Rcheck/` directory is `R CMD check` output; never commit edits to it.
- `.DS_Store` files appear under [R/](R/) and [vignettes/](vignettes/) — ignore, don't track.

## Open work

Tracked in [TODO.md](TODO.md).
