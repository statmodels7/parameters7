# Derivative Arrays of an Inverse Autoregressive Parameter

Assembles one order of \\\partial\Omega\\ for \\\Omega = U^\top
\mathrm{diag}(\tau) U\\ by the Leibniz rule, taken twice so that the
three factors cost \\2^m\\ products per component rather than \\3^m\\.

## Usage

``` r
ar_inv_derivative(s, eta, order)
```

## Arguments

- s:

  An
  [`AutoregressiveInvParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveInvParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A named list of `s@dimension` square matrices, keyed and ordered as
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
says.

## Details

Writing \\N = \mathrm{diag}(\tau)U\\, the inner Leibniz gives
\\\partial_T N = \sum\_{R \subseteq T} \mathrm{diag}(\partial_R \tau)
\partial\_{T \setminus R} U\\, and the outer one \\\partial_I \Omega =
\sum\_{S \subseteq I} (\partial_S U)^\top \partial\_{I \setminus S} N\\.
Both sums run over subsets of the index POSITIONS, which is what makes a
repeated index count with its multiplicity. Every \\\partial_S N\\ is
built once and read by every component that needs it.

## See also

[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
for the formula.
