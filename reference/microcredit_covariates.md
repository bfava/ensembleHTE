# Control Covariates for the Microcredit Dataset

A character vector of 24 covariate names used as controls in the
analysis of the Philippine microcredit experiment by Athey, Fava,
Karlan, Osman & Zinman (2025). Provided as a convenience for use in
formulas and subsetting.

## Usage

``` r
microcredit_covariates
```

## Format

A character vector of length 24 containing column names from the
[`microcredit`](https://bfava.com/ensembleHTE/reference/microcredit.md)
dataset: `"css_creditscorefinal"`, `"css_yearsataddress"`,
`"own_house"`, `"sari"`, `"own_anybus"`, `"max_yearsinbusiness"`,
`"css_regularworkers"`, `"css_traveltime"`, `"css_travelcost"`,
`"css_stockvalue"`, `"css_assetvalue"`, `"css_noofbusiness"`,
`"hhsize"`, `"age"`, `"gender"`, `"education"`, `"married"`,
`"savingsamt"`, `"hhinc_yrly_base"`, `"profits_yrly_base"`,
`"rev_yrly_base"`, `"exp_yrly_base"`, `"hhexp_yrly_base"`,
`"hhasset_yrly_base"`.

## Source

Karlan, D. and Zinman, J. (2011). Microcredit in Theory and Practice:
Using Randomized Credit Scoring for Impact Evaluation. *Science*,
**332**(6035), 1278–1284.
[doi:10.1126/science.1200138](https://doi.org/10.1126/science.1200138)

Athey, S., Fava, B., Karlan, D., Osman, A. and Zinman, J. (2025).
Profits and Social Impacts: Complements vs. Tradeoffs for Lenders in
Three Countries. Working paper.

## See also

[`microcredit`](https://bfava.com/ensembleHTE/reference/microcredit.md)

## Examples

``` r
data(microcredit_covariates)
microcredit_covariates
#>  [1] "css_creditscorefinal" "css_yearsataddress"   "own_house"           
#>  [4] "sari"                 "own_anybus"           "max_yearsinbusiness" 
#>  [7] "css_regularworkers"   "css_traveltime"       "css_travelcost"      
#> [10] "css_stockvalue"       "css_assetvalue"       "css_noofbusiness"    
#> [13] "hhsize"               "age"                  "gender"              
#> [16] "education"            "married"              "savingsamt"          
#> [19] "hhinc_yrly_base"      "profits_yrly_base"    "rev_yrly_base"       
#> [22] "exp_yrly_base"        "hhexp_yrly_base"      "hhasset_yrly_base"   

# Use to subset the microcredit data
data(microcredit)
dat <- microcredit[, c("hhinc_yrly_end", "treat", microcredit_covariates)]
```
