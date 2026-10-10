# Reject an Order That Would Difference a Difference

Signals an error when
[`param_d3logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.matrix_parameter.md)
or
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
is called for a family whose
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
is itself numerical. Both fallbacks difference
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md).
A family's own
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
is exact, and so is the base-class one **given**
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md);
in either case the differencing is one layer. Where the base-class
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
reads numerical arrays it is two, which is the nesting that the toolkit
avoids everywhere.

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

  The order requested, 3 or 4, which the message names.

## Value

Invisibly `TRUE`. For a family whose
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
comes from the base class and whose
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
or
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
comes from the base class too, an error that names the missing method
and the remedy.

## Details

With numerical arrays, the fourth order can come back larger than the
quantity that it estimates, and
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reports `TRUE` for the order whether the arrays are analytic or not, so
a consumer could not distinguish such a number from a usable one.

A family that writes its own
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
passes whether or not its
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
are analytic.

## See also

[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md),
which reports the routes read here, and
[`param_d3logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.matrix_parameter.md)
and
[`param_d4logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md),
the two callers.
