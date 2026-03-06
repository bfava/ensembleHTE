# Microcredit Experiment in the Philippines

Data from a randomized controlled trial that evaluated the impact of
individual-liability microloans on borrowers in the Philippines. Loans
were randomly assigned to applicants depending on credit score. The
dataset contains 1113 observations and 51 variables, prepared following
the processing in Athey, Fava, Karlan, Osman & Zinman (2025).

## Usage

``` r
microcredit
```

## Format

A data frame with 1113 rows and 51 variables:

- agg_id:

  Unique identifier for each individual.

- treat:

  Treatment indicator (1 = received microloan offer, 0 = control).

- prop_score:

  Propensity score (probability of treatment assignment).

- fu_survey_year:

  Year of the follow-up survey (2007 or 2008).

- hhinc_yrly_end:

  Household income at endline (yearly, Philippine pesos).

- profits_yrly_end:

  Business profits at endline (yearly, Philippine pesos).

- rev_yrly_end:

  Business revenue at endline (yearly, Philippine pesos).

- exp_yrly_end:

  Business expenses at endline (yearly, Philippine pesos).

- hhexp_yrly_end:

  Household expenses at endline (yearly, Philippine pesos).

- loan_size:

  Size of the microloan received (Philippine pesos; 0 if no loan taken).

- prinpaid:

  Principal amount paid back.

- interest_rate:

  Annual interest rate on the loan (percent).

- repaytime:

  Repayment time indicator.

- intpaid:

  Interest paid on the loan.

- penalty:

  Penalty fees paid.

- paid_total:

  Total amount paid (principal + interest + penalties).

- net_revenue:

  Net revenue earned by the lender from the loan.

- css_creditscorefinal:

  Credit score at application.

- css_yearsataddress:

  Years at current address at application.

- own_house:

  Indicator for home ownership (1 = yes).

- sari:

  Indicator for ownership of a sari-sari store (small convenience store;
  1 = yes).

- own_anybus:

  Indicator for owning any business (1 = yes).

- max_yearsinbusiness:

  Maximum years in business among owned businesses.

- css_regularworkers:

  Number of regular workers at application.

- css_traveltime:

  Travel time to lender branch (minutes).

- css_travelcost:

  Travel cost to lender branch (Philippine pesos).

- css_stockvalue:

  Stock value of the business at application (Philippine pesos).

- css_assetvalue:

  Asset value of the business at application (Philippine pesos).

- css_noofbusiness:

  Number of businesses at application.

- lower_window:

  Indicator for being in the lower window of the credit score
  distribution (1 = yes).

- app_year:

  Year of the loan application.

- hhsize:

  Household size.

- age:

  Age of the applicant (years).

- gender:

  Gender of the applicant (1 = female, 0 = male).

- education:

  Education level (categorical).

- married:

  Marital status (1 = married, 0 = not married).

- savingsamt:

  Savings amount (Philippine pesos).

- hhinc_yrly_base:

  Household income at baseline (yearly, Philippine pesos).

- profits_yrly_base:

  Business profits at baseline (yearly, Philippine pesos).

- rev_yrly_base:

  Business revenue at baseline (yearly, Philippine pesos).

- exp_yrly_base:

  Business expenses at baseline (yearly, Philippine pesos).

- hhexp_yrly_base:

  Household expenses at baseline (yearly, Philippine pesos).

- hhasset_yrly_base:

  Household assets at baseline (yearly, Philippine pesos).

- revenue_fees:

  Revenue from fees earned by the lender.

- revenue_interest:

  Revenue from interest earned by the lender.

- log_hhinc:

  Log of household income at endline.

- default_value:

  Default value of the loan (amount unpaid).

- high_school:

  Indicator for high school completion (1 = yes).

- college:

  Indicator for college completion (1 = yes).

- bank_profits_level:

  Bank profits in levels (Philippine pesos).

- bank_profits_pp:

  Bank profits per peso lent.

## Source

Karlan, D. and Zinman, J. (2011). Microcredit in Theory and Practice:
Using Randomized Credit Scoring for Impact Evaluation. *Science*,
**332**(6035), 1278–1284.
[doi:10.1126/science.1200138](https://doi.org/10.1126/science.1200138)

Athey, S., Fava, B., Karlan, D., Osman, A. and Zinman, J. (2025).
Profits and Social Impacts: Complements vs. Tradeoffs for Lenders in
Three Countries. Working paper.

## Details

The original experiment randomly assigned individual-liability
microloans to applicants at a microlender in the Philippines using
credit scoring as the randomization device. Applicants near the
credit-score cutoff were randomly assigned to treatment (loan approved)
or control (loan denied). The data include baseline covariates
(demographics, business characteristics, credit score components),
treatment assignment, loan repayment outcomes, and follow-up survey
outcomes measured 11–22 months after treatment.

## Examples

``` r
data(microcredit)
head(microcredit)
#>   agg_id treat prop_score fu_survey_year hhinc_yrly_end profits_yrly_end
#> 1 100519     1  0.8434874           2007          15750            15750
#> 2 101793     1  0.5465839           2008          66990            66990
#> 3 102790     1  0.8434874           2007           9450             9450
#> 4 103223     1  0.8434874           2007          29025            29025
#> 5 104470     1  0.8434874           2007           3750             3750
#> 6 105906     1  0.8434874           2007              0                0
#>   rev_yrly_end exp_yrly_end hhexp_yrly_end loan_size prinpaid interest_rate
#> 1        21750        29475         3000.0         0       NA            NA
#> 2        78750        11760         4500.0         0       NA            NA
#> 3        11250         1800        11250.0       625    10000            30
#> 4        30000          975        10500.0         0       NA            NA
#> 5         7500         9000         5812.5         0       NA            NA
#> 6            0            0        11250.0       625    10000            30
#>   repaytime intpaid penalty paid_total net_revenue css_creditscorefinal
#> 1         1      NA      NA         NA     0.00000                   48
#> 2         1      NA      NA         NA     0.00000                   42
#> 3         1  758.29       0   10758.29    47.39313                   60
#> 4         1      NA      NA         NA     0.00000                   52
#> 5         1      NA      NA         NA     0.00000                   54
#> 6         1  758.29       0   10758.29    47.39313                   56
#>   css_yearsataddress own_house sari own_anybus max_yearsinbusiness
#> 1                 44         1    1          1                   5
#> 2                  4         1    0          1                   4
#> 3                 15         1    1          1                  15
#> 4                 39         1    1          1                  10
#> 5                 39         1    1          1                   4
#> 6                 15         1    0          1                   4
#>   css_regularworkers css_traveltime css_travelcost css_stockvalue
#> 1                  0             10             30          25000
#> 2                  0              5             12         150000
#> 3                  0             20             24          17000
#> 4                  0             10              0          13800
#> 5                  0             10             24          15000
#> 6                  4             20              0           5000
#>   css_assetvalue css_noofbusiness lower_window app_year hhsize age gender
#> 1         140800                3            0     2006      1  44      1
#> 2              0                1            1     2006      1  46      0
#> 3          24000                2            0     2006      5  36      1
#> 4          11000                2            0     2006      4  39      0
#> 5          25000                1            0     2006      3  39      1
#> 6          15000                3            0     2006      6  38      1
#>   education married savingsamt hhinc_yrly_base profits_yrly_base rev_yrly_base
#> 1         1       0          0       21443.750         128662.50      224018.8
#> 2         2       1          0       57031.250         342187.50      342187.5
#> 3         2       1          0       10950.000          65700.00      123187.5
#> 4         1       1          0       21907.604         131445.62      250549.7
#> 5         1       1          0        5703.125          34218.75       48362.5
#> 6         1       1          0       23572.917         141437.50      234284.4
#>   exp_yrly_base hhexp_yrly_base hhasset_yrly_base revenue_fees revenue_interest
#> 1      95356.25         3390.00            1562.5    0.0000000          0.00000
#> 2          0.00         2325.00            1000.0    0.0000000          0.00000
#> 3      57487.50         2137.50           14012.5    0.6465497         46.74658
#> 4     119104.06         1050.00            2125.0    0.0000000          0.00000
#> 5      14143.75         2756.25            2875.0    0.0000000          0.00000
#> 6      92846.88         8421.00            7875.0    0.6465497         46.74658
#>   log_hhinc default_value high_school college bank_profits_level
#> 1  9.664659             0           1       1            0.00000
#> 2 11.112314             0           1       0            0.00000
#> 3  9.153876             0           1       0          -16.93826
#> 4 10.275947             0           1       1            0.00000
#> 5  8.229778             0           1       1            0.00000
#> 6  0.000000             0           1       1          -16.93826
#>   bank_profits_pp
#> 1         0.00000
#> 2         0.00000
#> 3         1.06774
#> 4         0.00000
#> 5         0.00000
#> 6         1.06774
summary(microcredit$treat)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0000  1.0000  1.0000  0.8005  1.0000  1.0000 
```
