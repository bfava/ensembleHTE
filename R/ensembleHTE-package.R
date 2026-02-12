#' ensembleHTE: Ensemble Methods for Heterogeneous Treatment Effect Estimation
#'
#' @description
#' The ensembleHTE package provides tools for estimating heterogeneous treatment
#' effects (HTE) and general predictions using ensemble learning methods. It combines
#' multiple machine learning algorithms with repeated cross-fitting to improve
#' statistical power.
#'
#' @details
#' The package provides two main estimation functions:
#' \itemize{
#'   \item \code{\link{ensemble_hte}}: Estimates individual treatment effects (ITE)
#'     using metalearner strategies (R-learner, T-learner, S-learner, X-learner)
#'   \item \code{\link{ensemble_pred}}: Estimates outcome predictions using an
#'     ensemble of ML algorithms
#' }
#'
#' Key features:
#' \itemize{
#'   \item Ensemble learning combining multiple ML algorithms (grf, lm, ranger, etc.)
#'   \item Repeated K-fold cross-fitting (M repetitions) for improved power
#'   \item Cross-validated Best Linear Predictor (BLP) weights for ensemble combination
#'   \item Analysis tools: GATES, BLP, CLAN, GAVS for characterizing treatment effects
#' }
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{ensemble_hte}}: Fit ensemble HTE model with metalearners
#'   \item \code{\link{ensemble_pred}}: Fit ensemble prediction model
#'   \item \code{\link{gates}}: Group Average Treatment Effects
#'   \item \code{\link{blp}}: Best Linear Predictor of CATE
#'   \item \code{\link{clan}}: Classification Analysis by covariate
#'   \item \code{\link{gavs}}: Group Averages
#' }
#'
#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @import mlr3
#' @import mlr3learners
#' @importFrom rlang ensym as_string
#' @importFrom stats predict coef cor median pnorm qnorm quantile sd var
#' @importFrom utils modifyList head
#' @importFrom methods is
NULL

# Suppress R CMD check notes for data.table non-standard evaluation
# These are column names used in data.table's non-standard evaluation
utils::globalVariables(c(
  # data.table column names
  ".", "..controls", "..covariate_vars", "..variables",
  "W1", "W2", "ci_lower", "ci_upper", "coef", "diff_top_all",
  "diff_top_bottom", "diff_top_else", "estimate", "fold", "group",
  "mean_all", "mean_bottom", "mean_else", "mean_top", "p_value",
  "p_value_top_all", "p_value_top_bottom", "p_value_top_else",
  "se", "se_diff", "se_top_all", "se_top_bottom", "se_top_else",
  "stars", "stars_ta", "stars_tb", "stars_te", "strategy",
  "t_value", "term", "variable", "weight"
))

#' @importFrom stats cov pnorm qnorm rnorm rbinom model.matrix var sd median
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Suppress mlr3/bbotk INFO logging to avoid verbose output during training/tuning
  lgr::get_logger("mlr3")$set_threshold("warn")
  lgr::get_logger("bbotk")$set_threshold("warn")
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {

  packageStartupMessage(
    "ensembleHTE v", utils::packageVersion("ensembleHTE"), " (beta)\n",
    "This package is under active development.\n",
    "Please report bugs or suggestions: https://github.com/bfava/ensembleHTE/issues"
  )
}
