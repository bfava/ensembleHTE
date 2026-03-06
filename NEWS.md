# ensembleHTE 0.1.1

## New features

* `Y`, `X`, and `D` in `ensemble_hte()` (and `Y`/`X` in `ensemble_pred()`) now
  accept column name(s) in addition to raw data. Pass `Y` as a single string,
  `X` as a character vector, or `D` as a single string together with `data` to
  have columns extracted automatically (e.g.,
  `Y = "hhinc_yrly_end", X = microcredit_covariates, D = "treat", data = microcredit`).
  Similarly, `prop_score` in `ensemble_hte()` now accepts a column name string
  (e.g., `prop_score = "prop_score"`).
* Prettier startup message with clickable links (in RStudio) using the `cli`
  package. `ensemble_news()` and `citation("ensembleHTE")` are now clickable
  in the loading message.
* Added `cli` to Imports.

## Improvements

* Significance stars now shown in `summary.ensemble_hte_fit()` (BLP beta1/beta2),
  `summary.ensemble_pred_fit()` (BLP intercept/slope), and restricted comparison
  print methods (`print.gates_restricted_results()`, `print.gavs_restricted_results()`).
* Signif. codes legend in restricted comparison print methods updated to the
  standard R format (`0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1`).
* `print.clan_results()` now computes column widths dynamically from the data,
  ensuring headers, values, standard errors, and separator lines are always
  properly aligned regardless of value magnitude.
* ITE Distribution section removed from `summary.ensemble_hte_fit()` output.
  Individual ITE point estimates are not directly interpretable; the summary
  now focuses on BLP and GATES inference.
* Prediction Summary section removed from `summary.ensemble_pred_fit()` output
  for the same reason. The returned object no longer includes `prediction_summary`.
* Internal formula is now rebuilt to reflect the actual covariates used, improving
  consistency when `Y ~ .` with exclusions is passed.
* `as.data.table()` conversion is now skipped when the input is already a
  `data.table`, avoiding unnecessary copies.

## Bug fixes

* Fixed `clan()` data subsetting when column names differ from the default
  internal names (`Y`, `D`, etc.), which could occur with the new column-name
  interface.
* Improved error messages when `Y`, `X`, or `D` are accidentally `NULL` (e.g.,
  from passing a non-existent column like `df$wrong_name`). The error now
  identifies which argument is `NULL` and suggests checking the column name.

## Documentation

* README substantially rewritten with a streamlined Quick Start using real
  `microcredit` data and the column-name interface, including an
  `ensemble_pred()` example.
* New walkthrough vignette (`vignettes/articles/microcredit-walkthrough.Rmd`)
  demonstrating a full analysis with the microcredit dataset.
* `microcredit` dataset documentation updated with Karlan & Zinman (2011)
  citation and DOI.

# ensembleHTE 0.1.0

* Initial CRAN release.
* Core estimation functions:
  - `ensemble_hte()`: Ensemble heterogeneous treatment effect estimation with
    metalearners (R-learner, T-learner, S-learner, X-learner).
  - `ensemble_pred()`: Ensemble prediction model for general outcome prediction.
  - `combine_ensembles()`: Combine multiple ensemble fits from distributed
    computing.
* Analysis functions:
  - `gates()`: Group Average Treatment Effects (GATES) estimation.
  - `blp()`: Best Linear Predictor of CATE.
  - `clan()`: Classification Analysis by covariates.
  - `blp_pred()`: Best Linear Predictor for prediction tasks.
  - `gavs()`: Group Averages for prediction tasks.
* Comparison functions:
  - `gates_restricted()`: Compare GATES between unrestricted and restricted
    ranking.
  - `gavs_restricted()`: Compare GAVS between unrestricted and restricted
    ranking.
* S3 methods for `print`, `summary`, and `plot`.
* Comprehensive documentation and vignettes.
* Test suite.
