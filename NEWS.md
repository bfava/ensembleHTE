# ensembleHTE 0.3.1

## New features

* `ensemble_hte()` and `ensemble_pred()` gain a `balance_folds_by` argument for
  stratified fold assignment. In `ensemble_hte()` the folds are already
  stratified so the treated share is balanced across the K folds (and, when a
  training subset is supplied via `train_idx`, so that training and
  prediction-only observations are spread evenly across folds); `balance_folds_by`
  additionally guarantees that a categorical variable (e.g. community/village,
  gender) is balanced across the K folds in every one of the M repetitions.
  Because the balancing variables are crossed with treatment, the balance is
  *joint* — e.g. with `balance_folds_by = "female"` the treated share *within
  each gender* is balanced across folds, not merely the marginal shares. Accepts
  a quoted column name, a character vector of column names, a variable holding a
  name, or a length-n vector. A message reports the number of stratum cells, and
  a warning fires when some cell has fewer than `K` observations (so balance is
  only approximate). Fit objects expose `$balance_folds_by`, and `print()` shows
  the balancing variables. This is distinct from `individual_id` (keep a unit's
  rows in one fold) and `se_cluster_id` (SE clustering).

# ensembleHTE 0.3.0

## New features

* `ensemble_hte()` and `ensemble_pred()` gain a `se_cluster_id` argument that
  decouples the level of cluster-robust standard errors from the level used to
  form cross-fitting folds. `individual_id` now controls only how folds are
  built (all of a unit's rows stay in one fold), while `se_cluster_id` controls
  the clustering of standard errors in the downstream analyses (`blp()`,
  `gates()`, `clan()`, `gavs()`, ...). If only one of the two is supplied, it is
  used for both roles. When a fit starts, an informative message reports the
  fold-splitting level and the SE-clustering level. Fit objects now expose
  `$se_cluster_id` (and the `$individual_id_name` / `$se_cluster_id_name`
  labels), and `print()` shows the two levels separately. A typical use is
  predicting outcomes for unobserved units within observed clusters (split by
  individual, cluster SEs by village).

## Bug fixes

* Panel data (`individual_id`) is now respected in the cross-validated
  ensemble step, not just the outer cross-fitting split. Previously, when
  `ensemble_strategy = "cv"`, the inner ensemble-weight CV split ignored the
  cluster identifier, so observations from the same individual could land in
  different ensemble folds and appear on both the fitting and prediction sides
  of that inner CV. The ensemble split (and the baseline ensemble that reuses
  it) in `ensemble_hte()` and `ensemble_pred()` now cluster by
  `individual_id`, matching the outer split.

# ensembleHTE 0.2.1

## Bug fixes

* `blp()` and `blp_pred()` no longer error with "subscript out of bounds"
  when a fit produces a constant (degenerate) prediction. Previously, a
  constant fitted CATE made the heterogeneity regressor collinear, so `lm()`
  dropped it and the hard-coded coefficient lookup failed. The affected
  coefficient (heterogeneity for `blp()`, slope for `blp_pred()`) now returns
  `NA` with an informative warning explaining that no heterogeneity (or
  predictive signal) was detected, rather than failing silently.

# ensembleHTE 0.2.0

## New features

* `ensemble_hte()` gains a `store_baseline` argument that controls whether
  predicted baseline outcomes are stored in the returned object.
  `"ensemble"` (default) stores an ensembled E[Y(0)|X] (T/S/X-learner) or
  E[Y|X] (R-learner) as a data.table with one column per repetition;
  `"all"` stores the full per-algorithm prediction array (n x A x M);
  `"none"` skips storage.
* `blp()`, `gates()`, and `gates_restricted()` gain a `baseline_as_control`
  argument. When `NULL` (default), the stored baseline is automatically
  included as a regression control if available. Including the predicted
  baseline absorbs residual outcome variance not attributable to the
  treatment, reducing standard errors in the BLP/GATES regressions and
  improving statistical power for detecting treatment effect heterogeneity.
  Set to `FALSE` to omit even if stored, or `TRUE` to require it (errors
  if no baseline was stored in the fit).

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
