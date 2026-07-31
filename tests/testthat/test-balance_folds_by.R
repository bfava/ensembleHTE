test_that("create_folds spreads each stratum cell near-evenly across folds", {
  set.seed(42)
  n <- 300
  female <- rep(0:1, length.out = n)
  D <- rbinom(n, 1, 0.5)
  train_idx <- rep(TRUE, n)

  stratify <- interaction(D, train_idx, female, drop = TRUE)
  folds <- create_folds(n, M = 1, K = 3, stratify_var = stratify)[[1]]

  # Within each (D x train x female) cell, the K fold counts differ by <= 1,
  # i.e. treatment share within each gender is balanced across folds (joint).
  tab <- table(stratify, folds)
  spread <- apply(tab, 1, function(r) max(r) - min(r))
  expect_true(all(spread <= 1))
})


test_that("ensemble_hte records, reports, and prints balance_folds_by", {
  skip_on_cran()

  set.seed(1)
  n <- 200
  dat <- data.frame(
    female = rep(0:1, length.out = n),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  expect_message(
    fit <- ensemble_hte(
      Y ~ X1 + X2, treatment = "D", data = dat,
      balance_folds_by = "female",
      M = 2, K = 2, algorithms = "lm", metalearner = "t"
    ),
    "balancing cross-fitting folds jointly by treatment and 'female'"
  )

  expect_identical(fit$balance_folds_by, "female")

  out <- capture.output(print(fit))
  expect_true(any(grepl("Folds balanced by:\\s+female", out)))
})


test_that("balance_folds_by accepts multiple columns and balances jointly with treatment", {
  skip_on_cran()

  set.seed(7)
  n <- 400
  dat <- data.frame(
    female = rep(0:1, length.out = n),
    region = rep(1:2, length.out = n),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  fit <- suppressMessages(ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    balance_folds_by = c("female", "region"),
    M = 1, K = 2, algorithms = "lm", metalearner = "t"
  ))

  expect_identical(fit$balance_folds_by, c("female", "region"))

  # Treated share within each (female, region) cell is balanced across folds.
  fold <- fit$splits[[1]]
  cell <- interaction(dat$female, dat$region, dat$D, drop = TRUE)
  tab <- table(cell, fold)
  spread <- apply(tab, 1, function(r) max(r) - min(r))
  expect_true(all(spread <= 1))
})


test_that("balance_folds_by validates inputs", {
  skip_on_cran()

  set.seed(2)
  n <- 100
  dat <- data.frame(
    D = rbinom(n, 1, 0.5), X1 = rnorm(n), X2 = rnorm(n),
    female = rep(0:1, length.out = n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  expect_error(
    ensemble_hte(
      Y ~ X1 + X2, treatment = "D", data = dat,
      balance_folds_by = "nonexistent",
      M = 1, K = 2, algorithms = "lm", metalearner = "t"
    ),
    "not found in data"
  )

  dat$bad <- dat$female
  dat$bad[1] <- NA
  expect_error(
    ensemble_hte(
      Y ~ X1 + X2, treatment = "D", data = dat,
      balance_folds_by = "bad",
      M = 1, K = 2, algorithms = "lm", metalearner = "t"
    ),
    "contains NA"
  )
})


test_that("balance_folds_by works with ensemble_pred", {
  skip_on_cran()

  set.seed(3)
  n <- 200
  dat <- data.frame(
    female = rep(0:1, length.out = n),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + rnorm(n)

  expect_message(
    fit <- ensemble_pred(
      Y ~ X1 + X2, data = dat, balance_folds_by = "female",
      M = 1, K = 2, algorithms = "lm"
    ),
    "balancing cross-fitting folds"
  )
  expect_identical(fit$balance_folds_by, "female")
})


test_that("balance_folds_by warns when stratum cells are too small for K", {
  set.seed(9)
  n <- 60
  strat <- factor(rep(1:20, each = 3))
  D <- rbinom(n, 1, 0.5)
  train_idx <- rep(TRUE, n)
  sv <- interaction(D, train_idx, strat, drop = TRUE)

  expect_warning(
    suppressMessages(.report_balance_folds(sv, "strat", K = 5)),
    "cannot be spread across all K folds"
  )
})


# ---------------------------------------------------------------------------
# "Does it actually do what we want?" — marginal balance + backward compat
# ---------------------------------------------------------------------------

test_that("balance_folds_by balances the marginal share across folds", {
  skip_on_cran()

  set.seed(11)
  n <- 400
  dat <- data.frame(
    female = rep(0:1, length.out = n),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  fit <- suppressMessages(ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    balance_folds_by = "female",
    M = 1, K = 4, algorithms = "lm", metalearner = "t"
  ))

  fold <- fit$splits[[1]]
  # Each fold's female share is within a small tolerance of the overall share.
  shares <- tapply(dat$female, fold, mean)
  expect_true(max(abs(shares - mean(dat$female))) < 0.05)

  # Fold sizes themselves stay balanced.
  sizes <- as.integer(table(fold))
  expect_true(max(sizes) - min(sizes) <= 4)
})

test_that("omitting balance_folds_by is backward compatible", {
  skip_on_cran()

  set.seed(12)
  n <- 160
  dat <- data.frame(
    D = rbinom(n, 1, 0.5), X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  # No message about balancing folds should be emitted.
  expect_no_message(
    fit <- ensemble_hte(
      Y ~ X1 + X2, treatment = "D", data = dat,
      M = 1, K = 2, algorithms = "lm", metalearner = "t"
    ),
    message = "balancing cross-fitting folds"
  )
  expect_null(fit$balance_folds_by)

  # Treatment is still balanced across folds by default.
  fold <- fit$splits[[1]]
  tab <- table(dat$D, fold)
  spread <- apply(tab, 1, function(r) max(r) - min(r))
  expect_true(all(spread <= 1))

  # print() does not show a "Folds balanced by" line.
  out <- capture.output(print(fit))
  expect_false(any(grepl("Folds balanced by", out)))
})


# ---------------------------------------------------------------------------
# Interactions with other fold-related choices
# ---------------------------------------------------------------------------

test_that("balance_folds_by + individual_id keeps units together and balances at unit level", {
  skip_on_cran()

  set.seed(13)
  nv <- 40; per <- 4; n <- nv * per
  dat <- data.frame(
    id = rep(seq_len(nv), each = per),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  # Balancing variable constant within unit.
  dat$female <- rep(rep(0:1, length.out = nv), each = per)
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  fit <- suppressMessages(ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    individual_id = "id", balance_folds_by = "female",
    M = 1, K = 2, algorithms = "lm", metalearner = "t"
  ))

  fold <- fit$splits[[1]]
  # Every unit's rows land in exactly one fold (keep-together preserved).
  units_per_fold <- tapply(fold, dat$id, function(z) length(unique(z)))
  expect_true(all(units_per_fold == 1))

  # Female share balanced across folds at the unit level.
  first <- !duplicated(dat$id)
  shares <- tapply(dat$female[first], fold[first], mean)
  expect_true(max(abs(shares - mean(dat$female[first]))) < 0.1)

  expect_identical(fit$balance_folds_by, "female")
})

test_that("balance_folds_by + se_cluster_id stores all three roles and runs downstream", {
  skip_on_cran()

  set.seed(14)
  nv <- 30; per <- 5; n <- nv * per
  dat <- data.frame(
    village = rep(seq_len(nv), each = per),
    id = seq_len(n),
    D = rbinom(n, 1, 0.5),
    female = rep(0:1, length.out = n),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  fit <- suppressMessages(ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    individual_id = "id", se_cluster_id = "village",
    balance_folds_by = "female",
    M = 1, K = 2, algorithms = "lm", metalearner = "t"
  ))

  expect_identical(fit$individual_id, dat$id)
  expect_identical(fit$se_cluster_id, dat$village)
  expect_identical(fit$balance_folds_by, "female")

  # Downstream analyses (which use the SE-clustering id) still run.
  expect_no_error(suppressWarnings(gates(fit, n_groups = 2)))
  expect_no_error(suppressWarnings(blp(fit)))
})

test_that("balance_folds_by respects a training subset (train_idx)", {
  skip_on_cran()

  set.seed(15)
  n <- 240
  dat <- data.frame(
    female = rep(0:1, length.out = n),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)
  train_idx <- rep(c(TRUE, TRUE, TRUE, FALSE), length.out = n)

  fit <- suppressMessages(ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    train_idx = train_idx, balance_folds_by = "female",
    M = 1, K = 3, algorithms = "lm", metalearner = "t"
  ))

  fold <- fit$splits[[1]]
  # Each fold contains both training and holdout observations.
  tab <- table(train_idx, fold)
  expect_true(all(tab > 0))
  # ITEs are produced for all observations (train and holdout).
  expect_false(anyNA(fit$ite[[1]]))
})

test_that("balance_folds_by works with ensemble_strategy = 'average'", {
  skip_on_cran()

  set.seed(16)
  n <- 160
  dat <- data.frame(
    female = rep(0:1, length.out = n),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  fit <- suppressMessages(ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    balance_folds_by = "female", ensemble_strategy = "average",
    M = 1, K = 2, algorithms = c("lm", "grf"), metalearner = "t"
  ))
  expect_identical(fit$balance_folds_by, "female")
  expect_false(anyNA(fit$ite[[1]]))
})


# ---------------------------------------------------------------------------
# Input-form robustness
# ---------------------------------------------------------------------------

test_that("balance_folds_by accepts a raw length-n vector and a factor", {
  skip_on_cran()

  set.seed(17)
  n <- 160
  dat <- data.frame(
    region = factor(rep(c("a", "b", "c"), length.out = n)),
    D = rbinom(n, 1, 0.5),
    X1 = rnorm(n), X2 = rnorm(n)
  )
  dat$Y <- dat$X1 + dat$D + rnorm(n)

  # Raw vector (not a column name) -> label derived from the expression.
  female_vec <- rep(0:1, length.out = n)
  fit_vec <- suppressMessages(ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    balance_folds_by = female_vec,
    M = 1, K = 2, algorithms = "lm", metalearner = "t"
  ))
  expect_length(fit_vec$balance_folds_by, 1L)

  # Factor column balances across folds.
  fit_fac <- suppressMessages(ensemble_hte(
    Y ~ X1 + X2, treatment = "D", data = dat,
    balance_folds_by = "region",
    M = 1, K = 2, algorithms = "lm", metalearner = "t"
  ))
  fold <- fit_fac$splits[[1]]
  # Balance is guaranteed per full stratum cell (region x treatment), since the
  # balancing variable is crossed with treatment.
  cell <- interaction(dat$region, dat$D, drop = TRUE)
  tab <- table(cell, fold)
  spread <- apply(tab, 1, function(r) max(r) - min(r))
  expect_true(all(spread <= 1))
})

test_that("balance_folds_by works with the matrix (Y, X, D) interface", {
  skip_on_cran()

  set.seed(18)
  n <- 160
  X <- data.frame(X1 = rnorm(n), X2 = rnorm(n))
  D <- rbinom(n, 1, 0.5)
  female <- rep(0:1, length.out = n)
  Y <- X$X1 + D + rnorm(n)

  fit <- suppressMessages(ensemble_hte(
    Y = Y, X = X, D = D,
    balance_folds_by = female,
    M = 1, K = 2, algorithms = "lm", metalearner = "t"
  ))
  expect_length(fit$balance_folds_by, 1L)
  expect_false(anyNA(fit$ite[[1]]))
})
