## Script to prepare the `microcredit` dataset and `microcredit_covariates`
## Source: Karlan & Zinman (2011), processed following Athey, Fava, Karlan,
##         Osman & Zinman (2025).

microcredit <- read.csv("data-raw/microcredit.csv")

## Covariate selection from Athey, Fava, Karlan, Osman & Zinman (2025)
microcredit_covariates <- c(
  "css_creditscorefinal", "css_yearsataddress", "own_house",
  "sari", "own_anybus", "max_yearsinbusiness",
  "css_regularworkers", "css_traveltime", "css_travelcost",
  "css_stockvalue", "css_assetvalue", "css_noofbusiness",
  "hhsize", "age", "gender",
  "education", "married", "savingsamt",
  "hhinc_yrly_base", "profits_yrly_base", "rev_yrly_base",
  "exp_yrly_base", "hhexp_yrly_base", "hhasset_yrly_base"
)

usethis::use_data(microcredit, microcredit_covariates, overwrite = TRUE)
