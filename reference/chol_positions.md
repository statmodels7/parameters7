# Positions of the Free Values in the Lower Triangle

Returns the row and column of each free value of a log-Cholesky
parameter, in the order the free vector uses: the \\p\\ diagonal entries
first, then the strictly below-diagonal entries column by column.
Everything in the family reads it, so the ordering is decided once here
and nowhere else.

## Usage

``` r
chol_positions(p)
```

## Arguments

- p:

  The side of the matrix, a single positive integer. At `p = 1` there
  are no below-diagonal entries and all three vectors have length 1.

## Value

A list with three vectors of length \\p(p+1)/2\\: the integer `row` and
`col` of each free value, and the logical `on_diagonal`.

## Details

The ordering is part of the interface, since `free_names` follows it and
consumers build their parameter tables from those names. At \\p = 3\\
the rows are `1 2 3 2 3 3`, the columns `1 2 3 1 1 2`, and `on_diagonal`
is `TRUE TRUE TRUE FALSE FALSE FALSE`, which labels the free values
`log_L1 log_L2 log_L3 L2.1 L3.1 L3.2`.

## See also

[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
which turns this into `free_names`, and
[`chol_assemble()`](https://statmodels7.github.io/parameters7/reference/chol_assemble.md),
which uses it to place the values in \\L\\.
