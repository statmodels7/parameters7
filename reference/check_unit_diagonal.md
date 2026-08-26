# Refuse a Correlation Block That Carries a Scale

Reads the diagonal of `correlation`'s own value at two probe free
vectors and signals an error unless every entry is 1. This is the
property \\\Sigma = D R D\\ rests on: the standard deviations are read
off \\D\\, so a block with a scale of its own describes the same matrix
from a whole ray of free vectors and the composite has one free value
too many.

## Usage

``` r
check_unit_diagonal(correlation)
```

## Arguments

- correlation:

  The
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  handed to
  [`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md),
  already checked for its class, its side and its rank.

## Value

Invisibly `TRUE`. A block whose diagonal departs from 1 at either probe
throws, with the family named and the worst entry quoted. A probe the
family cannot evaluate is skipped rather than reported, the constructor
having no standing to decide what a foreign chart admits.

## Why the probes are not the zero vector

A scale-carrying family is almost always written on a log link, so at a
zero free vector its scale is \\\exp(0) = 1\\ and its diagonal is
exactly the one this looks for. Measured over the six shipped families
of side four, a zero probe passes all five that should be rejected –
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
and
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
– while any non-zero probe catches every one of them, the worst diagonal
entry departing from 1 by between 0.35 and 8.27. The probes are
therefore `0.3, 0.4, ...` and a constant `-0.4`, and they are fixed
rather than drawn, so a rejection is reproducible.

Two probes cannot prove a property that is quantified over the whole
free space, and this does not claim to. What it catches is a family that
carries a scale, which is the way the requirement is broken in practice;
a family contrived to have a unit diagonal at exactly these two points
passes.
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)'s
round trip is what reports the ray itself.

## See also

[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md),
the caller, and
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
the shipped family with the property.
