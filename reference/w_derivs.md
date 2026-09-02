# Derivatives of the Reciprocal of One Minus a Squared Correlation

Returns the value and four derivatives of \\\eta \mapsto w = (1 -
\rho(\eta)^2)^{-1}\\, the quantity every inverse autoregression is
written in: it is the reciprocal of the innovation variance's share at
one lag, and it multiplies every entry of an inverse AR(1) and every
diagonal of an inverse AR(\\q\\).

## Usage

``` r
w_derivs(v)
```

## Arguments

- v:

  A list of five numbers: \\\rho\\ and its four derivatives in the free
  value, as
  [`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md)
  returns for one link.

## Value

A list of five numbers: \\w\\ and its four derivatives in the free
value.

## Details

Written by partial fractions, \\w = \tfrac{1}{2}\\(1-\rho)^{-1} +
(1+\rho)^{-1}\\\\, every order is an exact expression rather than a
repeated quotient rule,

\$\$w^{(k)} = \frac{k!}{2}\left\\(1-\rho)^{-(k+1)} + (-1)^k
(1+\rho)^{-(k+1)}\right\\,\$\$

and it is finite throughout \\\|\rho\| \< 1\\, which the correlation's
link guarantees. The result is then chained onto the free value.

## See also

[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
and
[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md),
which share it.
