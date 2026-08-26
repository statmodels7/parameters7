# Derivative Components of a Correlation Parameter

Assembles a whole derivative order of \\R = LL^\top\\ by the Leibniz
rule, the factor's derivatives coming from
[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)
and the trigonometric tables of
[`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md).
The four methods
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
through
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
of this family are one call each to this function.

## Usage

``` r
corr_derivative(s, eta, order)
```

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A list of `choose(s@n_free + order - 1, order)` symmetric matrices keyed
as `param_tuple_names(s, order)` and in that order, each `s@dimension`
by `s@dimension` with a zero diagonal.

## Details

Every term of the Leibniz sum whose factor derivative is `NULL` is
skipped, and most are: a multiset spanning two rows of \\L\\ contributes
nothing. The diagonal of the result is zero at every order above zero,
the diagonal of \\R\\ being the constant 1.

## See also

[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)
for the factor's derivatives,
[`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md)
for the sum, and
[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
the same construction for the log-Cholesky family.
