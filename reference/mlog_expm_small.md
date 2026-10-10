# Exponential of a Small Upper Triangular Matrix

Scaling and squaring with a Taylor series, written out here for the tiny
bidiagonal matrices that
[`dd_exp()`](https://statmodels7.github.io/parameters7/reference/dd_exp.md)
builds. The matrices are at most 5 by 5, so a general-purpose matrix
exponential, and the dependency that it would bring, is not needed.

## Usage

``` r
mlog_expm_small(a)
```

## Arguments

- a:

  A square numeric matrix, upper triangular in every call the package
  makes, and at most 5 by 5.

## Value

Its exponential, a numeric matrix of the same shape.

## See also

[`dd_exp()`](https://statmodels7.github.io/parameters7/reference/dd_exp.md),
the only caller.
