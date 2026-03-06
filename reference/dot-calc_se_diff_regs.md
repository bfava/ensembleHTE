# Compute SE of difference between two regression coefficients

Computes the estimate and standard error of the difference between
corresponding coefficients from two regressions (reg1 - reg2).

## Usage

``` r
.calc_se_diff_regs(reg1, reg2, coef_name)
```

## Arguments

- reg1:

  Result from reg_from_formula (list with coef, model_matrix, residuals)

- reg2:

  Result from reg_from_formula (list with coef, model_matrix, residuals)

- coef_name:

  Name of the coefficient to compare

## Value

Named vector with 'estimate' and 'se'
