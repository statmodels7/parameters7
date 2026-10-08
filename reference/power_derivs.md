# Derivatives of a Power, for Composition

Returns the derivatives of orders 1 to `order` of \\r \mapsto r^{m}\\ at
\\r\\, ready to be the outer map of a
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md)
call. They are \\m(m-1)\cdots(m-k+1)\\r^{m-k}\\ and vanish beyond order
\\m\\, the power being a polynomial.
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md)
calls it once per distinct lag.

## Usage

``` r
power_derivs(r, m, order)
```

## Arguments

- r:

  The point, a single number or a vector.

- m:

  The exponent, a non-negative integer. At `m = 0` all are 0, the power
  being the constant 1; at `m = 2` the first four are `c(2r, 2, 0, 0)`.

- order:

  The highest order wanted, an integer from 0 to 4.

## Value

A list of `order` elements, the first to the `order`-th derivative, each
the shape of `r`.

## See also

[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md),
which chains these onto a link's derivatives, and
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md),
the caller.
