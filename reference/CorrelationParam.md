# Correlation Matrix Parameter

The S7 class of correlation matrices, symmetric positive definite with a
unit diagonal, in the spherical parametrization of Rapisarda, Brigo and
Mercurio (2007). The unit diagonal holds by construction, never by a
correction, so every free vector gives a genuine correlation matrix.

[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
builds one. The free values are angles carried onto the real line, and
`param_params` records the row and column each belongs to together with
the `bounded_link(0, pi)` that carries it.

## Usage

``` r
CorrelationParam(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list(),
  dimension = integer(0),
  rank = integer(0),
  null_basis = integer(0)
)
```

## Arguments

- param_name:

  A single character string naming the family.

- n_free:

  The length \\d\\ of the free vector: a single non-negative integer,
  agreeing with `length(free_names)`.

- free_names:

  A character vector of length `n_free`, one label per free value, in
  the order the free vector holds them. Must be unique.

- param_params:

  A list of whatever the family needs in order to evaluate itself, read
  only by that family's own methods.

- dimension:

  The side \\p\\ of the matrix: a single integer, no smaller than 1. The
  validator rejects a vector and a value below 1.

- rank:

  The rank of the matrix the family produces, a single integer in
  `0:dimension`. It is a property of the family, so a family whose value
  is positive definite at every \\\eta\\ declares \\p\\ here.

- null_basis:

  A `dimension` by `dimension - rank` numeric matrix whose columns are
  an orthonormal basis of the common null space. Use
  [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
  to obtain one, or `matrix(numeric(0), dimension, 0)` for a full-rank
  family. The validator rejects any other shape, and reports both the
  rank and the shape when the two disagree.

## Value

An object of class `CorrelationParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `row` and `col`,
the position each angle belongs to, and `link`, a
`linkfunctions7::bounded_link(lwr = 0, upr = pi)`. `n_free` is
\\p(p-1)/2\\ and `rank` is \\p\\.

## References

Rapisarda, F., Brigo, D. and Mercurio, F. (2007). Parameterizing
correlations: a geometric interpretation. *IMA Journal of Management
Mathematics* **18**, 55-73.

## See also

[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
the constructor,
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
to give this a diagonal scale and make it a covariance, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
s <- correlation_matrix(3)
S7::S7_inherits(s, CorrelationParam)
#> [1] TRUE

# The free values are angles: p(p-1)/2 of them, one per below-diagonal entry.
c(n_free = s@n_free, p_choose_2 = 3 * 2 / 2)
#>     n_free p_choose_2 
#>          3          3 
s@free_names
#> [1] "z2.1" "z3.1" "z3.2"

# Every free vector gives a unit diagonal exactly, not approximately.
diag(param_value(s, c(2.5, -3, 1.7)))
#> v1 v2 v3 
#>  1  1  1 
```
