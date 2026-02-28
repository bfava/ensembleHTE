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
