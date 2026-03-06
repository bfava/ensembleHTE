# ============================================================================
# Pre-compute ensemble fit objects for vignettes
# ============================================================================
#
# This script generates the pre-computed fit objects used by the package
# vignettes. Vignettes load these .rds files instead of running the
# (expensive) ML training during vignette build.
#
# Run this script whenever:
#   - The microcredit dataset changes
#   - The ensemble_hte() or ensemble_pred() API changes
#   - You want to refresh the vignette output
#
# Usage:
#   source("data-raw/precompute_vignette_fits.R")
# ============================================================================

library(ensembleHTE)
data(microcredit)

# ---------------------------------------------------------------------------
# Shared variables (must match the vignette code)
# ---------------------------------------------------------------------------
hte_covars <- c("css_creditscorefinal", "lower_window", "own_anybus",
                "max_yearsinbusiness", "css_assetvalue")

profit_covars <- c("css_creditscorefinal", "own_anybus", "css_assetvalue",
                   "age", "gender", "hhinc_yrly_base")

has_loan <- microcredit$treat == 1 & microcredit$loan_size > 0

# Create income quintiles (needed for gates_restricted / gavs_restricted)
microcredit$hhinc_quintile <- cut(
  microcredit$hhinc_yrly_base,
  breaks = quantile(microcredit$hhinc_yrly_base, probs = seq(0, 1, by = 0.2)),
  labels = paste0("Q", 1:5),
  include.lowest = TRUE
)

out_dir <- "vignettes/precomputed"

# ---------------------------------------------------------------------------
# Complete Guide fits
# ---------------------------------------------------------------------------

cat("Fitting: Complete Guide — main HTE model ...\n")
fit_hte <- ensemble_hte(
  Y = "exp_yrly_end",
  D = "treat",
  X = hte_covars,
  data       = microcredit,
  prop_score = "prop_score",
  algorithms = c("lm", "grf"),
  M = 2, K = 3
)
saveRDS(fit_hte, file.path(out_dir, "fit_hte.rds"))

cat("Fitting: Complete Guide — prediction model ...\n")
fit_pred <- ensemble_pred(
  Y    = "bank_profits_pp",
  X    = profit_covars,
  data = microcredit,
  train_idx  = has_loan,
  algorithms = c("lm", "grf"),
  M = 2, K = 3
)
saveRDS(fit_pred, file.path(out_dir, "fit_pred.rds"))

cat("Fitting: Complete Guide — subset HTE model ...\n")
fit_subset_hte <- ensemble_hte(
  Y = "exp_yrly_end",
  D = "treat",
  X = hte_covars,
  data       = microcredit,
  prop_score = "prop_score",
  train_idx  = microcredit$lower_window == 1,
  algorithms = c("lm", "grf"),
  M = 2, K = 3
)
saveRDS(fit_subset_hte, file.path(out_dir, "fit_subset_hte.rds"))

# ---------------------------------------------------------------------------
# Quick Guide fits (different covariate set = 4 covars, no lower_window)
# ---------------------------------------------------------------------------
hte_covars_quick <- c("css_creditscorefinal", "own_anybus",
                      "max_yearsinbusiness", "css_assetvalue")

cat("Fitting: Quick Guide — HTE model ...\n")
fit_hte_quick <- ensemble_hte(
  Y = "exp_yrly_end",
  D = "treat",
  X = hte_covars_quick,
  data       = microcredit,
  prop_score = "prop_score",
  algorithms = c("lm", "grf"),
  M = 2, K = 4
)
saveRDS(fit_hte_quick, file.path(out_dir, "fit_hte_quick.rds"))

cat("Fitting: Quick Guide — prediction model ...\n")
fit_pred_quick <- ensemble_pred(
  Y    = "bank_profits_pp",
  X    = c("css_creditscorefinal", "own_anybus", "css_assetvalue",
           "age", "gender", "hhinc_yrly_base"),
  data = microcredit,
  train_idx  = has_loan,
  algorithms = c("lm", "grf"),
  M = 2, K = 4
)
saveRDS(fit_pred_quick, file.path(out_dir, "fit_pred_quick.rds"))

cat("Done! Pre-computed fit files saved to:", out_dir, "\n")
cat("Files:\n")
cat(paste(" ", list.files(out_dir, pattern = "\\.rds$")), sep = "\n")
