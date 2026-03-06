# Parse Column Name from Various Input Types

Internal helper to parse a column name specification that can be
provided as:

- An unquoted symbol (e.g., `treatment = D`)

- A quoted string (e.g., `treatment = "D"`)

- A variable containing a string (e.g.,
  `my_var <- "D"; treatment = my_var`)

## Usage

``` r
parse_column_name(expr, env, arg_name = "variable", data = NULL)
```

## Arguments

- expr:

  The captured expression (from
  [`substitute()`](https://rdrr.io/r/base/substitute.html))

- env:

  The environment to evaluate in (typically
  [`parent.frame()`](https://rdrr.io/r/base/sys.parent.html))

- arg_name:

  Name of the argument (for error messages)

- data:

  Optional data.frame to validate column existence

## Value

Character string with the column name
