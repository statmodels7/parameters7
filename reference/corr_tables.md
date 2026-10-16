# Sines and Cosines of a Correlation Parameter's Angles

Returns, for every angle, the value and the derivatives to order `order`
in the free value of both \\\sin\theta\\ and \\\cos\theta\\. All the
derivatives of the family are built from these tables: an entry of \\L\\
is a product of such factors, so differentiating it replaces each factor
by the derivative of the matching order, and the tables contain every
derivative that is required.

## Usage

``` r
corr_tables(s, eta, order)
```

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object, whose `param_params$link` is read.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The highest derivative order to compute, an integer from 0 to 4.

## Value

A list with two components, `sin` and `cos`, each a list of `s@n_free`
lists of `order + 1` numbers: the value at index 1 and the derivatives
in the free value at the following indices.

## Details

Each angle depends on exactly one free value, so a table is indexed by
free value and cross terms do not appear at this level. The chain from
\\\mathrm{d}/\mathrm{d}\theta\\ to \\\mathrm{d}/\mathrm{d}\eta\\ is
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)'s,
applied to the link's own derivatives, so the accuracy is the link's and
the derivatives are computed without differencing. Order 0 gives only
the value, and the derivatives are not evaluated in that case.

## See also

[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)
for the chain rule,
[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)
which reads these tables, and
[`corr_logdet_chains()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_chains.md),
which reads the `sin` half.
