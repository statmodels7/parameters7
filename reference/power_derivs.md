# Derivatives of a Power, for Composition

Returns the four derivatives of \\r \mapsto r^{m}\\ at \\r\\, ready to
be the outer map of a
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
call. They are \\m(m-1)\cdots(m-k+1)\\r^{m-k}\\ and vanish beyond order
\\m\\, the power being a polynomial.
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md)
calls it once per distinct lag.

## Usage

``` r
power_derivs(r, m)
```

## Arguments

- r:

  The point, a single number or a vector.

- m:

  The exponent, a non-negative integer. At `m = 0` all four are 0, the
  power being the constant 1; at `m = 2` they are `c(2r, 2, 0, 0)`.

## Value

A list of four elements, the first to fourth derivative, each the shape
of `r`.

## See also

[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md),
which chains these onto a link's derivatives, and
[`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md),
the caller.
