---
description: "Use when editing R source files: roxygen2 docblocks, exports, internal helpers, S3 methods, package imports. Enforces ensembleHTE conventions for documentation, naming, and regeneration."
applyTo: "R/**/*.R"
---
# roxygen2 + R Package Conventions (ensembleHTE)

## Documentation

- Every function gets a roxygen2 block. Order tags: `@title` (implicit first line) → `@description` → `@param` → `@return` → `@examples` → `@export` or `@keywords internal` → `@references`.
- Exported functions need `@export` AND at least one runnable `@examples` block. Wrap slow examples in `\donttest{}`; never use `\dontrun{}` for normal usage.
- Internal helpers MUST have `@keywords internal` (and no `@export`) so they don't pollute the package index.
- Use `@inheritParams` to avoid duplicating common params (`M`, `K`, `data`, `algorithms`, etc.) across `ensemble_hte()` / `ensemble_pred()` / analysis functions.
- Math: use `\eqn{...}` for inline, `\deqn{...}` for display. Cite Fava (2025) via `@references` — see [R/ensemble_hte.R](../../R/ensemble_hte.R) for the canonical pattern.

## Naming & file layout

- Internal helpers are prefixed with `.` (e.g. `.fit_model_grf`, `.blp_single`, `.encode_factors_if_needed`) and live in the same file as their exported caller — do not create a new `R/` file for them.
- S3 methods: name `method.class` (e.g. `print.ensemble_hte_fit`, `plot.gates_results`) and tag with `@export` — roxygen registers `S3method(...)` in [NAMESPACE](../../NAMESPACE) automatically.
- New result classes: assign `class(x) <- c("<name>_results", "list")` and provide at minimum `print.*` and ideally `plot.*`.

## Imports

- Never use `library()` or `require()` inside `R/`. Use `@importFrom pkg fn` for ≤3 symbols from a package; use `pkg::fn()` inline for one-offs.
- Hard dependencies live in `DESCRIPTION` `Imports:`; `glmnet`, `gbm`, `ranger`, `xgboost` stay in `Suggests:` and must be guarded with `requireNamespace("pkg", quietly = TRUE)` before use.
- Algorithm dispatch: `"grf"` → `grf::causal_forest` (see [R/ml_training.R](../../R/ml_training.R)); everything else routes through `mlr3learners`. Don't add a new branch without updating both `.fit_model_*` and `.predict_model`.

## Validation & messaging

- Use `cli::cli_abort()` / `cli::cli_warn()` / `cli::cli_inform()` instead of `stop()` / `warning()` / `message()` for user-facing errors. Plain `stop()` is acceptable only inside `.`-prefixed helpers for invariants.
- Validate user input at the top of exported functions; let internal helpers assume valid input.

## After editing

1. Run `devtools::document()` to regenerate `man/*.Rd` and `NAMESPACE`. Never hand-edit either.
2. Run `devtools::load_all(".")` then `devtools::test()` for affected modules.
3. Before committing public-API changes, regenerate vignette fits: `source("data-raw/precompute_vignette_fits.R")` (and `precompute_walkthrough.R` for the complete guide).

## Don't

- Don't add `@param`/`@return` for arguments you didn't touch — leave existing prose alone.
- Don't write to `man/` or `NAMESPACE` by hand.
- Don't introduce `magrittr` `%>%`; package already uses base R and `data.table` — match the surrounding file.
