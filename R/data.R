#' Microcredit Experiment in the Philippines
#'
#' Data from a randomized controlled trial that evaluated the impact of
#' individual-liability microloans on borrowers in the Philippines. Loans
#' were randomly assigned to applicants depending on credit score. The dataset 
#' contains 1113 observations and 51 variables, prepared following the 
#' processing in Athey, Fava, Karlan, Osman & Zinman (2025).
#'
#' @format A data frame with 1113 rows and 51 variables:
#' \describe{
#'   \item{agg_id}{Unique identifier for each individual.}
#'   \item{treat}{Treatment indicator (1 = received microloan offer, 0 =
#'     control).}
#'   \item{prop_score}{Propensity score (probability of treatment assignment).}
#'   \item{fu_survey_year}{Year of the follow-up survey (2007 or 2008).}
#'   \item{hhinc_yrly_end}{Household income at endline (yearly, Philippine
#'     pesos).}
#'   \item{profits_yrly_end}{Business profits at endline (yearly, Philippine
#'     pesos).}
#'   \item{rev_yrly_end}{Business revenue at endline (yearly, Philippine
#'     pesos).}
#'   \item{exp_yrly_end}{Business expenses at endline (yearly, Philippine
#'     pesos).}
#'   \item{hhexp_yrly_end}{Household expenses at endline (yearly, Philippine
#'     pesos).}
#'   \item{loan_size}{Size of the microloan received (Philippine pesos; 0 if no
#'     loan taken).}
#'   \item{prinpaid}{Principal amount paid back.}
#'   \item{interest_rate}{Annual interest rate on the loan (percent).}
#'   \item{repaytime}{Repayment time indicator.}
#'   \item{intpaid}{Interest paid on the loan.}
#'   \item{penalty}{Penalty fees paid.}
#'   \item{paid_total}{Total amount paid (principal + interest + penalties).}
#'   \item{net_revenue}{Net revenue earned by the lender from the loan.}
#'   \item{css_creditscorefinal}{Credit score at application.}
#'   \item{css_yearsataddress}{Years at current address at application.}
#'   \item{own_house}{Indicator for home ownership (1 = yes).}
#'   \item{sari}{Indicator for ownership of a sari-sari store (small
#'     convenience store; 1 = yes).}
#'   \item{own_anybus}{Indicator for owning any business (1 = yes).}
#'   \item{max_yearsinbusiness}{Maximum years in business among owned
#'     businesses.}
#'   \item{css_regularworkers}{Number of regular workers at application.}
#'   \item{css_traveltime}{Travel time to lender branch (minutes).}
#'   \item{css_travelcost}{Travel cost to lender branch (Philippine pesos).}
#'   \item{css_stockvalue}{Stock value of the business at application
#'     (Philippine pesos).}
#'   \item{css_assetvalue}{Asset value of the business at application
#'     (Philippine pesos).}
#'   \item{css_noofbusiness}{Number of businesses at application.}
#'   \item{lower_window}{Indicator for being in the lower window of the credit
#'     score distribution (1 = yes).}
#'   \item{app_year}{Year of the loan application.}
#'   \item{hhsize}{Household size.}
#'   \item{age}{Age of the applicant (years).}
#'   \item{gender}{Gender of the applicant (1 = female, 0 = male).}
#'   \item{education}{Education level (categorical).}
#'   \item{married}{Marital status (1 = married, 0 = not married).}
#'   \item{savingsamt}{Savings amount (Philippine pesos).}
#'   \item{hhinc_yrly_base}{Household income at baseline (yearly, Philippine
#'     pesos).}
#'   \item{profits_yrly_base}{Business profits at baseline (yearly, Philippine
#'     pesos).}
#'   \item{rev_yrly_base}{Business revenue at baseline (yearly, Philippine
#'     pesos).}
#'   \item{exp_yrly_base}{Business expenses at baseline (yearly, Philippine
#'     pesos).}
#'   \item{hhexp_yrly_base}{Household expenses at baseline (yearly, Philippine
#'     pesos).}
#'   \item{hhasset_yrly_base}{Household assets at baseline (yearly, Philippine
#'     pesos).}
#'   \item{revenue_fees}{Revenue from fees earned by the lender.}
#'   \item{revenue_interest}{Revenue from interest earned by the lender.}
#'   \item{log_hhinc}{Log of household income at endline.}
#'   \item{default_value}{Default value of the loan (amount unpaid).}
#'   \item{high_school}{Indicator for high school completion (1 = yes).}
#'   \item{college}{Indicator for college completion (1 = yes).}
#'   \item{bank_profits_level}{Bank profits in levels (Philippine pesos).}
#'   \item{bank_profits_pp}{Bank profits per peso lent.}
#' }
#'
#' @details
#' The original experiment randomly assigned individual-liability microloans to
#' applicants at a microlender in the Philippines using credit scoring as the
#' randomization device. Applicants near the credit-score cutoff were randomly
#' assigned to treatment (loan approved) or control (loan denied). The data
#' include baseline covariates (demographics, business characteristics, credit
#' score components), treatment assignment, loan repayment outcomes, and
#' follow-up survey outcomes measured 11--22 months after treatment.
#'
#' @source
#' Karlan, D. and Zinman, J. (2011). Microcredit in Theory and Practice: Using
#' Randomized Credit Scoring for Impact Evaluation. \emph{Science},
#' \strong{332}(6035), 1278--1284.
#' \doi{10.1126/science.1200138}
#'
#' Athey, S., Fava, B., Karlan, D., Osman, A. and Zinman, J. (2025). Profits
#' and Social Impacts: Complements vs. Tradeoffs for Lenders in Three
#' Countries. Working paper.
#'
#' @examples
#' data(microcredit)
#' head(microcredit)
#' summary(microcredit$treat)
"microcredit"

#' Control Covariates for the Microcredit Dataset
#'
#' A character vector of 24 covariate names used as controls in the analysis
#' of the Philippine microcredit experiment by Athey, Fava, Karlan, Osman &
#' Zinman (2025). Provided as a convenience for use in formulas and subsetting.
#'
#' @format A character vector of length 24 containing column names from the
#'   \code{\link{microcredit}} dataset:
#'   \code{"css_creditscorefinal"}, \code{"css_yearsataddress"},
#'   \code{"own_house"}, \code{"sari"}, \code{"own_anybus"},
#'   \code{"max_yearsinbusiness"}, \code{"css_regularworkers"},
#'   \code{"css_traveltime"}, \code{"css_travelcost"},
#'   \code{"css_stockvalue"}, \code{"css_assetvalue"},
#'   \code{"css_noofbusiness"}, \code{"hhsize"}, \code{"age"},
#'   \code{"gender"}, \code{"education"}, \code{"married"},
#'   \code{"savingsamt"}, \code{"hhinc_yrly_base"},
#'   \code{"profits_yrly_base"}, \code{"rev_yrly_base"},
#'   \code{"exp_yrly_base"}, \code{"hhexp_yrly_base"},
#'   \code{"hhasset_yrly_base"}.
#'
#' @source
#' Karlan, D. and Zinman, J. (2011). Microcredit in Theory and Practice: Using
#' Randomized Credit Scoring for Impact Evaluation. \emph{Science},
#' \strong{332}(6035), 1278--1284.
#' \doi{10.1126/science.1200138}
#'
#' Athey, S., Fava, B., Karlan, D., Osman, A. and Zinman, J. (2025). Profits
#' and Social Impacts: Complements vs. Tradeoffs for Lenders in Three
#' Countries. Working paper.
#'
#' @seealso \code{\link{microcredit}}
#'
#' @examples
#' data(microcredit_covariates)
#' microcredit_covariates
#'
#' # Use to subset the microcredit data
#' data(microcredit)
#' dat <- microcredit[, c("hhinc_yrly_end", "treat", microcredit_covariates)]
"microcredit_covariates"
