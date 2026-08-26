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

The scale term is \\p \log \gamma_0\\, whose derivatives in \\\gamma_0\\
are \\(-1)^{k-1}(k-1)!/\gamma_0^{k}\\, and each partial autocorrelation
contributes \\(p-k)\log(1 - r_k^2)\\, whose derivatives come from
[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md)
applied to the two factors \\1 - r\\ and \\1 + r\\. Both are then
carried onto the free scale by
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md),
the Faa di Bruno chain with a one-dimensional inner map.

## See also

[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
and
[`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md)
for the two pieces, and
[`param_dlogdet.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md),
which calls this.
