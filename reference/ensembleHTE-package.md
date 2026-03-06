# ensembleHTE: Ensemble Methods for Heterogeneous Treatment Effect Estimation

The ensembleHTE package provides tools for estimating heterogeneous
treatment effects (HTE) and general predictions using ensemble learning
methods. It combines multiple machine learning algorithms with repeated
cross-fitting to improve statistical power.

## Details

The package provides two main estimation functions:

- [`ensemble_hte`](https://bfava.github.io/ensembleHTE/reference/ensemble_hte.md):
  Estimates individual treatment effects (ITE) using metalearner
  strategies (R-learner, T-learner, S-learner, X-learner)

- [`ensemble_pred`](https://bfava.github.io/ensembleHTE/reference/ensemble_pred.md):
  Estimates outcome predictions using an ensemble of ML algorithms

Key features:

- Ensemble learning combining multiple ML algorithms (grf, lm, ranger,
  etc.)

- Repeated K-fold cross-fitting (M repetitions) for improved power

- Cross-validated Best Linear Predictor (BLP) weights for ensemble
  combination

- Analysis tools: GATES, BLP, CLAN, GAVS for characterizing treatment
  effects

## Main functions

- [`ensemble_hte`](https://bfava.github.io/ensembleHTE/reference/ensemble_hte.md):
  Fit ensemble HTE model with metalearners

- [`ensemble_pred`](https://bfava.github.io/ensembleHTE/reference/ensemble_pred.md):
  Fit ensemble prediction model

- [`gates`](https://bfava.github.io/ensembleHTE/reference/gates.md):
  Group Average Treatment Effects

- [`blp`](https://bfava.github.io/ensembleHTE/reference/blp.md): Best
  Linear Predictor of CATE

- [`clan`](https://bfava.github.io/ensembleHTE/reference/clan.md):
  Classification Analysis by covariate

- [`gavs`](https://bfava.github.io/ensembleHTE/reference/gavs.md): Group
  Averages

## See also

Useful links:

- <https://bfava.github.io/ensembleHTE/>

- <https://github.com/bfava/ensembleHTE>

- Report bugs at <https://github.com/bfava/ensembleHTE/issues>

## Author

**Maintainer**: Bruno Fava <brunovnfava@gmail.com>
([ORCID](https://orcid.org/0009-0002-8218-516X))
