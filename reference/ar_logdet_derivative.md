# Log-Determinant Components of an Autoregressive Parameter

Assembles one derivative order of the log-determinant. It is a sum with
one term per free value, so a mixed component is exactly zero and is
returned without arithmetic; only the \\q + 1\\ pure components are
computed.

## Usage

``` r
ar_logdet_derivative(s, eta, order)
```

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A numeric vector of `choose(s@n_free + order - 1, order)` values keyed
as `param_tuple_names(s, order)` and in that order.

## Details

The scale term is \\p \log \gamma_0\\: the derivatives of \\\log
\gamma_0\\ in \\\gamma_0\\ are \\(-1)^{k-1}(k-1)!/\gamma_0^{k}\\, and
the factor \\p\\ is applied afterwards; they are carried onto the free
scale by
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md),
the Faa di Bruno chain with a one-dimensional inner map. Each partial
autocorrelation contributes \\(p-k)\log(1 - r_k^2) =
(p-k)\log\mathrm{sech}^2 z_k\\ under the rhobit link, whose derivatives
in \\z_k\\ are \\p - k\\ times those of
[`log_sech2_deriv()`](https://statmodels7.github.io/parameters7/reference/sech2.md),
written in the free value with no chain.

## See also

[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)
and
[`log_sech2_deriv()`](https://statmodels7.github.io/parameters7/reference/sech2.md)
for the two pieces, and
[`param_dlogdet.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md),
which calls this.
