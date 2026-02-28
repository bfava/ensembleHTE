# ==============================================================================
# Tests for analysis functions: blp, blp_pred, clan, gates_restricted, gavs_restricted
# ==============================================================================

# Shared test data creation helper
create_test_data <- function(n = 200, seed = 123) {
  set.seed(seed)
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  X3 <- rbinom(n, 1, 0.5)  # Binary variable for restrict_by
  # Create a discrete version of X1 for restrict_by comparisons (3 groups)
  X1_group <- cut(X1, breaks = quantile(X1, probs = c(0, 1/3, 2/3, 1)), 
                  labels = c("low", "mid", "high"), include.lowest = TRUE)
  D <- rbinom(n, 1, 0.5)
  tau <- 1 + 0.5 * X1
  Y <- 0.5 * X1 + 0.3 * X2 + D * tau + rnorm(n, sd = 0.5)
  
  data.frame(Y = Y, D = D, X1 = X1, X2 = X2, X3 = X3, X1_group = X1_group)
}

# ==============================================================================
# Tests for blp()
# ==============================================================================

test_that("blp works with ensemble_hte_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  # Use T-learner which produces non-zero ITEs more reliably
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  result <- blp(fit)
  
  expect_s3_class(result, "blp_results")
  expect_true(!is.null(result$estimates))
  expect_true(!is.null(result$targeted_outcome))
  expect_equal(result$fit_type, "hte")
  
  # Print should work
  expect_output(print(result), "BLP")
})


test_that("blp works with ensemble_pred_fit and treatment", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  # Must specify treatment for pred fits
  result <- blp(fit, treatment = "D")
  
  expect_s3_class(result, "blp_results")
  expect_equal(result$fit_type, "pred")
  
  # Print should work
  expect_output(print(result), "BLP")
})


test_that("blp with ensemble_pred_fit requires treatment", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  expect_error(
    blp(fit),
    "treatment must be specified"
  )
})


test_that("blp accepts treatment as vector", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  result <- blp(fit, treatment = data$D)
  
  expect_s3_class(result, "blp_results")
})


# ==============================================================================
# Tests for blp_pred()
# ==============================================================================

test_that("blp_pred works with ensemble_pred_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  result <- blp_pred(fit)
  
  expect_s3_class(result, "blp_pred_results")
  expect_true(!is.null(result$estimates))
  expect_equal(result$fit_type, "pred")
  
  # Print should work
  expect_output(print(result), "BLP")
})


test_that("blp_pred works with ensemble_hte_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  result <- blp_pred(fit)
  
  expect_s3_class(result, "blp_pred_results")
  expect_equal(result$fit_type, "hte")
})


test_that("blp_pred subset options work with ensemble_hte_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Default (NULL) uses all observations
  result_all <- blp_pred(fit)
  expect_s3_class(result_all, "blp_pred_results")
  
  # Subset with logical vector
  result_subset <- blp_pred(fit, subset = rep(c(TRUE, FALSE), 75))
  expect_s3_class(result_subset, "blp_pred_results")
})


# ==============================================================================
# Tests for clan()
# ==============================================================================

test_that("clan works with ensemble_pred_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm")
  )
  
  # Specify variables to analyze
  result <- clan(fit, variables = c("X1", "X2"), n_groups = 3)
  
  expect_s3_class(result, "clan_results")
  expect_true(!is.null(result$estimates))
  expect_equal(result$n_groups, 3)
  
  # Print should work
  expect_output(print(result), "CLAN")
})


test_that("clan works with ensemble_hte_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  result <- clan(fit, variables = c("X1", "X2"), n_groups = 3)
  
  expect_s3_class(result, "clan_results")
})


test_that("clan plot method works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm")
  )
  
  result <- clan(fit, variables = c("X1"), n_groups = 3)
  
  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# Tests for gates_restricted()
# ==============================================================================

test_that("gates_restricted rejects continuous restrict_by variable", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Continuous X1 should be rejected with informative error
  expect_error(
    gates_restricted(fit, restrict_by = "X1", n_groups = 3),
    "continuous variable"
  )
})


test_that("gates_restricted works with ensemble_hte_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Use restrict_by parameter with categorical variable
  result <- gates_restricted(fit, restrict_by = "X1_group", n_groups = 3)
  
  expect_s3_class(result, "gates_restricted_results")
  expect_true(!is.null(result$unrestricted))
  expect_true(!is.null(result$restricted))
  
  # Print should work
  expect_output(print(result), "GATES")
})


test_that("gates_restricted works with ensemble_pred_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm")
  )
  
  result <- gates_restricted(fit, restrict_by = "X1_group", treatment = "D", n_groups = 3)
  
  expect_s3_class(result, "gates_restricted_results")
  expect_equal(result$fit_type, "pred")
})


test_that("gates_restricted with ensemble_pred_fit requires treatment", {
  skip_on_cran()
  
  data <- create_test_data(n = 100)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 2,
    algorithms = c("lm")
  )
  
  expect_error(
    gates_restricted(fit, restrict_by = "X1_group", n_groups = 3),
    "treatment must be specified"
  )
})


test_that("gates_restricted plot method works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  result <- gates_restricted(fit, restrict_by = "X1_group", n_groups = 3)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# Tests for gavs_restricted()
# ==============================================================================

test_that("gavs_restricted works with ensemble_pred_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm")
  )
  
  result <- gavs_restricted(fit, restrict_by = "X1_group", n_groups = 3)

  expect_s3_class(result, "gavs_restricted_results")
  expect_true(!is.null(result$unrestricted))
  expect_true(!is.null(result$restricted))

  # Print should work
  expect_output(print(result), "GAVS")
})


test_that("gavs_restricted works with ensemble_hte_fit", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
result <- gavs_restricted(fit, restrict_by = "X1_group", n_groups = 3)
  
  expect_s3_class(result, "gavs_restricted_results")
})


test_that("gavs_restricted plot method works", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm")
  )
  
  result <- gavs_restricted(fit, restrict_by = "X1_group", n_groups = 3)

  p <- plot(result)
  expect_s3_class(p, "ggplot")
})


# ==============================================================================
# Tests for subset argument in analysis functions
# ==============================================================================

test_that("gavs subset argument works with logical vector", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Use first 100 observations
  subset_logical <- c(rep(TRUE, 100), rep(FALSE, 50))
  
  result <- gavs(fit, n_groups = 3, subset = subset_logical)
  
  expect_s3_class(result, "gavs_results")
  expect_equal(result$n_used, 100)
})


test_that("gavs subset argument works with integer indices", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Use first 75 observations
  subset_int <- 1:75
  
  result <- gavs(fit, n_groups = 3, subset = subset_int)
  
  expect_s3_class(result, "gavs_results")
  expect_equal(result$n_used, 75)
})


test_that("gates subset argument works", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Use first 100 observations
  subset_int <- 1:100
  
  result <- gates(fit, n_groups = 3, subset = subset_int)
  
  expect_s3_class(result, "gates_results")
})


test_that("blp subset argument works", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Use first 100 observations
  subset_int <- 1:100
  
  result <- blp(fit, subset = subset_int)
  
  expect_s3_class(result, "blp_results")
})


test_that("clan subset argument works", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Use first 100 observations
  subset_logical <- c(rep(TRUE, 100), rep(FALSE, 50))
  
  result <- clan(fit, variables = c("X1", "X2"), n_groups = 3, subset = subset_logical)
  
  expect_s3_class(result, "clan_results")
})


test_that("clan group_on argument works", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  subset_logical <- c(rep(TRUE, 100), rep(FALSE, 50))
  
  # group_on = "all": groups from all obs, means on subset
  result_all <- clan(fit, variables = c("X1", "X2"), n_groups = 3,
                     subset = subset_logical, group_on = "all")
  expect_s3_class(result_all, "clan_results")
  expect_equal(result_all$group_on, "all")
  
  # group_on = "analysis": groups formed within the analysis subset
  result_sub <- clan(fit, variables = c("X1", "X2"), n_groups = 3,
                     subset = subset_logical, group_on = "analysis")
  expect_s3_class(result_sub, "clan_results")
  expect_equal(result_sub$group_on, "analysis")
  
  # The two modes should give different results when subsetting
  expect_false(isTRUE(all.equal(result_all$estimates$diff_top_bottom,
                                result_sub$estimates$diff_top_bottom)))
})


test_that("clan group_on has no effect without subsetting", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Without subset, both group_on modes should be identical
  result_all <- clan(fit, variables = c("X1", "X2"), n_groups = 3, group_on = "all")
  result_sub <- clan(fit, variables = c("X1", "X2"), n_groups = 3, group_on = "analysis")
  
  expect_equal(result_all$estimates$diff_top_bottom,
               result_sub$estimates$diff_top_bottom)
  expect_equal(result_all$estimates$mean_top,
               result_sub$estimates$mean_top)
})


test_that("clan group_on validates input", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Invalid group_on value should error
  expect_error(
    clan(fit, variables = c("X1", "X2"), group_on = "invalid"),
    "arg"
  )
})


test_that("subset argument validates input correctly", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # Wrong length logical vector should error
  expect_error(
    gavs(fit, n_groups = 3, subset = c(TRUE, FALSE)),
    "length"
  )
  
  # Out of bounds integer indices should error
  expect_error(
    gavs(fit, n_groups = 3, subset = 1:200),
    "indices"
  )
  
  # Invalid subset type should error
  expect_error(
    gavs(fit, n_groups = 3, subset = "invalid"),
    "must be"
  )
})


test_that("subset with ensemble_pred_fit works", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_pred(
    Y ~ X1 + X2,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm")
  )
  
  # Subset with integer indices
  subset_int <- 1:100
  
  result <- gavs(fit, n_groups = 3, subset = subset_int)
  expect_s3_class(result, "gavs_results")
  expect_equal(result$n_used, 100)
  
  # Subset with logical vector
  subset_logical <- rep(c(TRUE, FALSE), 75)
  result2 <- gavs(fit, n_groups = 3, subset = subset_logical)
  expect_s3_class(result2, "gavs_results")
  expect_equal(result2$n_used, 75)
})


test_that("subset='all' works correctly", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # subset="all" should use all observations
  result <- gavs(fit, n_groups = 3, subset = "all")
  
  expect_s3_class(result, "gavs_results")
  expect_equal(result$n_used, 150)
})


test_that("compare functions work with subset argument", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  subset_int <- 1:120
  
  # gavs_restricted with subset
  result_gavs <- gavs_restricted(fit, restrict_by = "X1_group", n_groups = 3, subset = subset_int)
  expect_s3_class(result_gavs, "gavs_restricted_results")

  # gates_restricted with subset
  result_gates <- gates_restricted(fit, restrict_by = "X1_group", n_groups = 3, subset = subset_int)
  expect_s3_class(result_gates, "gates_restricted_results")
})


# ---- GATES group_on tests ----

test_that("gates group_on argument works", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  subset_logical <- c(rep(TRUE, 100), rep(FALSE, 50))
  
  # group_on = "all": groups from all obs, regress on subset
  result_all <- gates(fit, n_groups = 3, subset = subset_logical, group_on = "all")
  expect_s3_class(result_all, "gates_results")
  expect_equal(result_all$group_on, "all")
  
  # group_on = "analysis": groups formed within the analysis subset
  result_sub <- gates(fit, n_groups = 3, subset = subset_logical, group_on = "analysis")
  expect_s3_class(result_sub, "gates_results")
  expect_equal(result_sub$group_on, "analysis")
  
  # The two modes should give different results when subsetting
  expect_false(isTRUE(all.equal(result_all$estimates$estimate,
                                result_sub$estimates$estimate)))
})


test_that("gates group_on has no effect without subsetting", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  result_all <- gates(fit, n_groups = 3, group_on = "all")
  result_sub <- gates(fit, n_groups = 3, group_on = "analysis")
  
  expect_equal(result_all$estimates$estimate, result_sub$estimates$estimate)
  expect_equal(result_all$top_bottom$estimate, result_sub$top_bottom$estimate)
})


test_that("gates group_on validates input", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  expect_error(
    gates(fit, n_groups = 3, group_on = "invalid"),
    "arg"
  )
})


# ---- GAVS group_on tests ----

test_that("gavs group_on argument works", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  subset_logical <- c(rep(TRUE, 100), rep(FALSE, 50))
  
  # group_on = "all": groups from all obs, regress on subset
  result_all <- gavs(fit, n_groups = 3, subset = subset_logical, group_on = "all")
  expect_s3_class(result_all, "gavs_results")
  expect_equal(result_all$group_on, "all")
  
  # group_on = "analysis": groups formed within the analysis subset
  result_sub <- gavs(fit, n_groups = 3, subset = subset_logical, group_on = "analysis")
  expect_s3_class(result_sub, "gavs_results")
  expect_equal(result_sub$group_on, "analysis")
  
  # The two modes should give different results when subsetting
  expect_false(isTRUE(all.equal(result_all$estimates$estimate,
                                result_sub$estimates$estimate)))
})


test_that("gavs group_on has no effect without subsetting", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  result_all <- gavs(fit, n_groups = 3, group_on = "all")
  result_sub <- gavs(fit, n_groups = 3, group_on = "analysis")
  
  expect_equal(result_all$estimates$estimate, result_sub$estimates$estimate)
  expect_equal(result_all$top_bottom$estimate, result_sub$top_bottom$estimate)
})


test_that("gavs group_on validates input", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  expect_error(
    gavs(fit, n_groups = 3, group_on = "invalid"),
    "arg"
  )
})


# ---- Small cell warnings ----

test_that("gates warns when fold x restrict_by cells are too small", {
  skip_on_cran()
  
  data <- create_test_data(n = 60)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # With n=60, K=3 folds (~20 per fold), and 5 restrict_by levels,
  # some fold x restrict_by cells will have ~4 obs but n_groups=5
  data$fine_strata <- as.factor(rep(1:5, length.out = 60))

  expect_warning(
    gates(fit, n_groups = 5, restrict_by = data$fine_strata),
    "fold.*restrict_by.*too small for"
  )
})

test_that("gavs warns when fold x restrict_by cells are too small", {
  skip_on_cran()
  
  data <- create_test_data(n = 60)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  data$fine_strata <- as.factor(rep(1:5, length.out = 60))
  
  expect_warning(
    gavs(fit, n_groups = 5, restrict_by = data$fine_strata),
    "fold.*restrict_by.*too small for"
  )
})

test_that("clan warns when fold cells are too small", {
  skip_on_cran()
  
  # Very small dataset with many folds -> some folds have <3 obs
  data <- create_test_data(n = 30)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # n=30, K=3 -> ~10 per fold; n_groups=15 -> 10 < 15, should warn
  expect_warning(
    clan(fit, variables = c("X1", "X2"), n_groups = 15),
    "fold.*too small for"
  )
})

test_that("no small cell warning when cells are adequate", {
  skip_on_cran()
  
  data <- create_test_data(n = 150)
  
  fit <- ensemble_hte(
    Y ~ X1 + X2,
    treatment = D,
    data = data,
    M = 2,
    K = 3,
    algorithms = c("lm"),
    metalearner = "t"
  )
  
  # n=150, K=3 -> ~50 per fold; n_groups=3 -> 50 >> 3, no warning
  expect_no_warning(gates(fit, n_groups = 3))
  expect_no_warning(gavs(fit, n_groups = 3))
  expect_no_warning(clan(fit, variables = c("X1", "X2"), n_groups = 3))
})
