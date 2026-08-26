# One Row of a Diagnostic Table

Builds one row of the table
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
accumulates and returns. Strings are kept as strings, so the assembled
table's `check` and `status` columns come out character.

## Usage

``` r
check_row(name, status, statistic = NA_real_)
```

## Arguments

- name:

  The name of the check, as it appears in the printed table:
  `"membership"`, `"round trip"`, `"first derivatives"` and the rest.

- status:

  One of `"OK"`, `"FAIL"` or `"NOT CHECKED"`. The third is used where
  the quantity comes from a numerical fallback and no independent route
  exists to compare against.

- statistic:

  The number the verdict rests on, a relative discrepancy, or `NA_real_`
  for a check that is structural and has no number, such as the shape
  and name check.

## Value

A one-row data frame with columns `check`, `status` and `statistic`.

## See also

[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
which assembles these rows.
