# Refuse an Order That Would Difference a Difference

Signals an error when the caller asks
[`param_d3logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.matrix_parameter.md)
or
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
of a family whose derivative arrays are themselves numerical. Both
fallbacks difference
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
which is an exact identity **given**
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md);
where those are supplied the differencing is one layer, and where they
are not it is two, which is the nesting the toolkit forbids everywhere.

## Usage

``` r
check_analytic_arrays(s, order)
```

## Arguments

- s:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  object.

- order:

  The order being asked for, 3 or 4, which the message names.

## Value

Invisibly `TRUE`. A family whose
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
or
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
comes from the base class throws instead, with both the missing method
and the remedy named.

## Details

The refusal is not a matter of accuracy alone. Measured on a \\4 \times
4\\ AR(1) covariance whose family supplies
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and nothing else, order three came back \\7.5 \times 10^{-3}\\ against a
quantity of size 4.45 and order four came back **9.07** against a
quantity of size 2.17, which is four times the size of what it
estimates. Nothing downstream could tell such a number from a usable
one:
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reports `TRUE` for the order in both regimes.

Only the two arrays are read, and not
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
itself. A family that writes its own second-order log-determinant out is
asked for nothing further, that method being reached before this guard
is.

## See also

[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md),
which answers the question this asks, and
[`param_d3logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.matrix_parameter.md)
and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md),
the two callers.
