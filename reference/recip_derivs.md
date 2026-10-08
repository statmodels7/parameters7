# Derivatives of a Reciprocal, for Composition

Returns the value and the derivatives of \\\eta \mapsto 1/h(\eta)\\ from
the value and the derivatives of \\h\\, by composing \\x^{-1}\\ onto
them, to the order `v` carries.

## Usage

``` r
recip_derivs(v)
```

## Arguments

- v:

  A list: the value of \\h\\ and its derivatives to some order up to
  four, as
  [`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md)
  returns for one link.

## Value

A list of the same length, for \\1/h\\.

## Details

[`power_derivs()`](https://statmodels7.github.io/parameters7/reference/power_derivs.md)
cannot serve here: it truncates beyond the exponent, which is right for
a non-negative integer power and wrong for \\-1\\, where every order is
non-zero. The outer derivatives written out are \\-x^{-2}\\,
\\2x^{-3}\\, \\-6x^{-4}\\ and \\24x^{-5}\\.

## See also

[`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md),
the caller, and
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)
for the chain.
