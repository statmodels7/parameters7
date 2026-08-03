# Positions of the Free Values in the Lower Triangle

The row and column of each free value of a log-Cholesky structure, in
the order the free vector uses: the diagonal first, then the
below-diagonal entries column by column.

## Usage

``` r
chol_positions(p)
```

## Arguments

- p:

  The side of the matrix.

## Value

A list with the integer vectors `row`, `col` and the logical
`on_diagonal`.

## Details

The ordering is fixed and is part of the contract, because `free_names`
depends on it and every consumer builds parameter tables from those
names.
