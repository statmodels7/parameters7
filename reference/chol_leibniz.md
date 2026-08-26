# Derivative Components of a Log-Cholesky Parameter

Assembles a whole derivative order of \\M = L L^\top\\ by the Leibniz
rule, in compiled code. Each component of \\\partial^S M\\ distributes
the differentiations of \\S\\ over the two factors, and since every
surviving \\\partial^T L\\ is a single-entry matrix (see
[`chol_dfactor()`](https://statmodels7.github.io/parameters7/reference/chol_dfactor.md)),
each term of the sum is one row, one column or one cell of a product.
The kernel walks those entries directly and never multiplies two
matrices.

The four methods
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
of this family are one call each to this function.

## Usage

``` r
chol_leibniz(s, eta, order)
```

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A list of `choose(s@n_free + order - 1, order)` symmetric matrices keyed
as `param_tuple_names(s, order)` and in that order, each `s@dimension`
by `s@dimension` with dimnames `v1`, `v2`, ...

## Details

Measured against
[`.chol_leibniz_r()`](https://statmodels7.github.io/parameters7/reference/dot-chol_leibniz_r.md),
the dense R twin, on the same free vector: **30x** at \\p = 8\\ and
order 4 (0.295 s against 8.89 s over 82251 components), 30x at \\p = 8\\
order 2, 35x at \\p = 5\\ order 4 and 19x at \\p = 3\\ order 2. The two
agree **exactly**, to 0, at \\p = 5\\ order 4.

## See also

[`.chol_leibniz_r()`](https://statmodels7.github.io/parameters7/reference/dot-chol_leibniz_r.md),
the R twin the tests hold it against,
[`chol_dfactor()`](https://statmodels7.github.io/parameters7/reference/chol_dfactor.md)
for the factor's derivatives, and
[`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md).
