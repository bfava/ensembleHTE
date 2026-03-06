# Get Significance Stars

Converts p-values to significance stars following standard convention:

- `***`: p \< 0.001

- `**`: p \< 0.01

- `*`: p \< 0.05

- `.`: p \< 0.1

- (empty): p \>= 0.1

## Usage

``` r
get_stars(p)
```

## Arguments

- p:

  Numeric vector of p-values

## Value

Character vector of significance stars

## Examples

``` r
if (FALSE) { # \dontrun{
get_stars(c(0.0001, 0.005, 0.03, 0.08, 0.5))
# Returns: "***" "**" "*" "." ""
} # }
```
