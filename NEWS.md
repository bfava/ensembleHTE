# Changelog

All notable changes to ensembleHTE will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-12-08

### Added
- Initial public release
- Core estimation functions:
  - `ensemble_hte()`: Ensemble heterogeneous treatment effect estimation with metalearners (R-learner, T-learner, S-learner, X-learner)
  - `ensemble_pred()`: Ensemble prediction model for general outcome prediction
  - `combine_ensembles()`: Combine multiple ensemble fits from distributed computing
- Analysis functions:
  - `gates()`: Group Average Treatment Effects (GATES) estimation
  - `blp()`: Best Linear Predictor of CATE
  - `clan()`: Classification Analysis by covariates
  - `blp_pred()`: Best Linear Predictor for prediction tasks
  - `gavs()`: Group Averages for prediction tasks
- Comparison functions:
  - `gates_compare()`: Compare GATES between unrestricted and restricted ranking
  - `gavs_compare()`: Compare GAVS between unrestricted and restricted ranking
- S3 methods for `print`, `summary`, and `plot`
- Comprehensive documentation and vignettes
- Test suite

[0.1.0]: https://github.com/bfava/ensembleHTE-dev/releases/tag/v0.1.0
