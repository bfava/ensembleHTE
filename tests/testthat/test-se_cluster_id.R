test_that("se_cluster_id decouples fold splitting from SE clustering", {
  skip_on_cran()

  set.seed(123)
  nv <- 25
  per <- 6
  n <- nv * per
  dat <- data.frame(
    village = rep(1:nv, each = per),
    id = 1:n,
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D * (1 + 0.5 * dat$X1) + rnorm(n, sd = 0.5)

  # Split by individual, cluster SEs by village
  expect_message(
    fit <- ensemble_hte(
      Y ~ X1 + X2, treatment = "D", data = dat,
      individual_id = "id", se_cluster_id = "village",
      M = 2, K = 2, algorithms = "lm", metalearner = "t"
    ),
    "splitting cross-fitting folds by 'id'.*clustering standard errors by 'village'"
  )

  # The two identifiers are stored separately with their labels
  expect_identical(fit$individual_id, dat$id)
  expect_identical(fit$se_cluster_id, dat$village)
  expect_identical(fit$individual_id_name, "id")
  expect_identical(fit$se_cluster_id_name, "village")

  # Downstream analysis uses the SE-clustering identifier and runs
  expect_no_error(gates(fit, n_groups = 2))
  expect_no_error(blp(fit))
})

test_that("supplying only one identifier uses it for both roles", {
  skip_on_cran()

  set.seed(1)
  nv <- 20
  per <- 5
  n <- nv * per
  dat <- data.frame(
    village = rep(1:nv, each = per),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n, sd = 0.5)

  # Only se_cluster_id
  fit_se <- ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    se_cluster_id = "village",
    M = 2, K = 2, algorithms = "lm", metalearner = "t"
  )
  expect_identical(fit_se$individual_id, dat$village)
  expect_identical(fit_se$se_cluster_id, dat$village)
  expect_identical(fit_se$individual_id_name, "village")
  expect_identical(fit_se$se_cluster_id_name, "village")

  # Only individual_id (backward-compatible behavior)
  fit_ind <- ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    individual_id = "village",
    M = 2, K = 2, algorithms = "lm", metalearner = "t"
  )
  expect_identical(fit_ind$se_cluster_id, dat$village)
})

test_that("ensemble_pred supports se_cluster_id", {
  skip_on_cran()

  set.seed(7)
  nv <- 20
  per <- 6
  n <- nv * per
  dat <- data.frame(
    village = rep(1:nv, each = per),
    id = 1:n,
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + rnorm(n, sd = 0.5)

  fp <- ensemble_pred(
    Y ~ X1 + X2, data = dat,
    individual_id = "id", se_cluster_id = "village",
    M = 2, K = 2, algorithms = "lm"
  )
  expect_identical(fp$individual_id, dat$id)
  expect_identical(fp$se_cluster_id, dat$village)
  expect_no_error(blp_pred(fp))
})
