#!/usr/bin/env Rscript
library(ensembleHTE)
data(microcredit)

hte_covars <- c("css_creditscorefinal", "lower_window", "own_anybus",
                "max_yearsinbusiness", "css_assetvalue")

set.seed(2026)
fit_hte_wt <- ensemble_hte(
  Y = "exp_yrly_end",
  D = "treat",
  X = hte_covars,
  data       = microcredit,
  prop_score = "prop_score",
  algorithms = c("lm", "grf"),
  M = 10,
  K = 4
)
saveRDS(fit_hte_wt, "vignettes/precomputed/fit_hte_walkthrough.rds")
cat("Saved fit_hte_walkthrough.rds\n")

profit_covars <- c("css_creditscorefinal", "lower_window", "own_anybus",
                   "max_yearsinbusiness", "css_assetvalue", "css_stockvalue",
                   "age", "gender", "hhinc_yrly_base", "exp_yrly_base")
has_loan <- microcredit$treat == 1 & microcredit$loan_size > 0

set.seed(2026)
fit_pred_wt <- ensemble_pred(
  Y    = "bank_profits_pp",
  X    = profit_covars,
  data = microcredit,
  train_idx  = has_loan,
  algorithms = c("lm", "grf"),
  M = 10,
  K = 4
)
saveRDS(fit_pred_wt, "vignettes/precomputed/fit_pred_walkthrough.rds")
cat("Saved fit_pred_walkthrough.rds\n")
