# Reject a Link That Does Not Carry the Whole Free Line to Positive Entries

Checks the two properties every family taking a link argument needs, and
signals an error naming the bounds otherwise: that the link's declared
range lies in the non-negative half line, and that it is defined on the
whole real line. Called at construction by
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md).

## Usage

``` r
check_positive_link(link)
```

## Arguments

- link:

  The object to check.

## Value

Invisibly `TRUE`. An object that is not a linkfunctions7 link throws
`'link' must be a linkfunctions7 link object.`; a link whose lower bound
is negative, or which is defined on part of the line, throws a message
quoting the offending bounds.

## Onto the positive half line

A diagonal entry of a positive definite matrix is positive, so a link
onto the whole real line would let a caller build a matrix outside the
set without anything saying so. The test is on the **lower** bound
alone, `b[1] >= 0`, so a link with a bounded range is accepted:
`logit_link()`, whose range is \\(0, 1)\\, produces a positive definite
matrix with entries below 1, which is a legitimate thing to want.

## And from the whole of it

The free vector is unconstrained by design, and every family's page
states that any vector in \\\mathbb{R}^d\\ gives a valid matrix. A link
defined on part of the line breaks that: `sqrt_link()`,
`inverse_link()`, `inverse_sq_link()` and `power_link()` at a positive
exponent all reach the positive entries from the positive predictors
alone, and
[`linkfunctions7::linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/linkinv.html)
of them is even, so `-2` and `2` give the same entry and the round trip
returns \\\lvert \eta \rvert\\.

Measured before the second test existed, with the default free-value
draws: `autoregressive(5, 1)`, `ar1(5)` and `compound_symmetry(4)` under
a square root link died inside
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
with `missing value where TRUE/FALSE needed`, and `diagonal_matrix(3)`
and `scalar_matrix(3)` with `system is computationally singular`,
neither message naming the link. Both properties are read at
construction, which is the only place they can be: at evaluation time no
free value is inadmissible.

## See also

[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
for the fullest statement of the two conditions, and
[`linkfunctions7::eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.html),
which answers the second.
