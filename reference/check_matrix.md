# Validate a Matrix Handed Back to a Structure

Checks that `m` is a square symmetric numeric matrix of the structure's
dimension.

## Usage

``` r
check_matrix(s, m, tol = 1e-08)
```

## Arguments

- s:

  A
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  object.

- m:

  The matrix supplied by the caller.

- tol:

  The relative tolerance for the symmetry check.

## Value

`m`, symmetrised.
