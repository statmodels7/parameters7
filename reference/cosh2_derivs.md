# Derivatives of the Reciprocal of One Minus a Squared Correlation

Returns the value and the derivatives to order `order` of \\w = (1 -
\rho^2)^{-1}\\ in the free value \\z\\, for \\\rho = \tanh z\\ under the
rhobit link. \\w\\ is the quantity every inverse autoregression is
written in: it is the reciprocal of the innovation variance's share at
one lag, and it multiplies every entry of an inverse AR(1) and every
diagonal of an inverse AR(\\q\\).

## Usage

``` r
cosh2_derivs(z, order)
```

## Arguments

- z:

  The free value, a single number.

- order:

  The highest order wanted, an integer from 0 to 5.

## Value

A list of `order + 1` numbers: \\w\\ and its derivatives in \\z\\.

## Details

Under the rhobit link \\w = \cosh^2 z = (1 + \cosh 2z)/2\\, so every
order is elementary: \\w^{(k)} = 2^{k-1}\sinh 2z\\ for odd \\k\\ and
\\2^{k-1}\cosh 2z\\ for even \\k\\. Written in \\z\\, no difference
\\1 - \rho^2\\ is formed, and the result keeps its accuracy where
\\\rho\\ rounds to \\-1\\ or 1.

## See also

[`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md),
whose factors call it, and
[`ar1_inv_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_inv_pattern.md),
which builds the inverse AR(1) pattern from it.
