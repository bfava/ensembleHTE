# Fetch update status from GitHub

Internal helper that fetches `inst/update_status.dcf` from the public
GitHub repository.

## Usage

``` r
.fetch_update_status()
```

## Value

A matrix from [`read.dcf()`](https://rdrr.io/r/base/dcf.html), or `NULL`
on failure.
