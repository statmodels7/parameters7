# Reject a Correlation Block That Carries a Scale

Reads the diagonal of `correlation`'s own value at two probe free
vectors and signals an error unless every entry is 1. This is the
property on which \\\Sigma = D R D\\ rests: the standard deviations are
read off \\D\\, so a block with a scale of its own describes the same
matrix from a whole ray of free vectors and the composite has one free
value too many.

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

Invisibly `TRUE`. For a block whose diagonal departs from 1 at either
probe, an error that names the family and quotes the worst entry. A
probe at which the family cannot be evaluated is skipped.

## The probe free vectors

A scale-carrying family is almost always written on a log link, so at a
zero free vector its scale is \\\exp(0) = 1\\ and its diagonal is 1. A
zero probe would therefore accept the scale-carrying families of the
package
([`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)),
while the non-zero probes reject each of them. The probes are
`0.3, 0.4, ...` and a constant `-0.4`, fixed and not drawn, so a
rejection is reproducible.

Two probes cannot prove a property that holds over the whole free space.
The check catches a family that carries a scale, which is the usual way
in which the requirement is broken; a family constructed to have a unit
diagonal at exactly these two points passes.
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)'s
round trip is what reports the ray itself.

## See also

[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md),
the caller, and
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
the shipped family with the property.
