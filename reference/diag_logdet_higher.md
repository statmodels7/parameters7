# Higher Derivatives of a Diagonal Log-Determinant

Assembles the derivative components of orders two to four of
\\\log\|M\|\\ for a diagonal family. The log-determinant is \\\sum_i
\log h(\eta_i)\\, a sum of functions of one free value each, so every
mixed component is exactly zero and a pure one is the matching
derivative from
[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md),
counted once per entry the free value owns.

## Usage

``` r
diag_logdet_higher(s, eta, order)
```

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 2, 3 or 4.

## Value

A numeric vector of `choose(s@n_free + order - 1, order)` entries, keyed
as `param_tuple_names(s, order)` and in that order. Under the log link
it is all zeros at every order above 1, \\\log h(\eta) = \eta\\ being
linear.

## See also

[`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
for the derivatives of \\\log h\\, and
[`param_d3logdet.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.DiagMatrixParam.md)
and
[`param_d4logdet.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.DiagMatrixParam.md),
the two callers.
