# The R Twin of the Compiled Leibniz Assembly

Computes the same components as `chol_leibniz_cpp`, through the dense
matrix products of
[`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md).
It is the independent reference against which the tests check that
kernel, so a change to one side that is not made to both shows up as a
disagreement. It is not called on a production path;
[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md)
is.

## Usage

``` r
.chol_leibniz_r(s, eta, order)
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
as `param_tuple_names(s, order)`, identical to
[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md)'s.

## Details

The two routes agree exactly, both summing the same Leibniz terms, with
the compiled one skipping only additions of structural zeros.

The compiled route is the production one because every derivative of the
factor is a single-entry matrix, so a Leibniz term is one row, one
column or one cell, where this twin forms a full \\p \times p\\ product.

## See also

[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
the compiled route this mirrors, and
[`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md),
the dense Leibniz sum it uses.
