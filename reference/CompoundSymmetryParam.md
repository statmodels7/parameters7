# Compound Symmetry Parameter

The S7 class of compound-symmetric covariance matrices: equal variances
and one common correlation, so **two** free values at every dimension.
It is the covariance of an exchangeable set of measurements, and the one
a random intercept induces.

[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
builds one. `param_params` records the two links, and the correlation's
is a `bounded_link(-1/(p-1), 1)` whose lower end depends on `dimension`:
two objects of different sizes carry different links.

## Usage

``` r
CompoundSymmetryParam(
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

An object of class `CompoundSymmetryParam`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `link_scale` and
`link_rho`. `n_free` is 2 at every `dimension`, and `rank` is \\p\\.

## See also

[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
the constructor,
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
for the other two-value family, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# Two free values whatever the dimension.
s <- compound_symmetry(3)
S7::S7_inherits(s, CompoundSymmetryParam)
#> [1] TRUE
vapply(2:6, function(p) compound_symmetry(p)@n_free, integer(1))
#> [1] 2 2 2 2 2

# The correlation's link depends on the dimension, its lower bound being
# -1/(p-1): below that the matrix would not be positive definite.
t(vapply(c(2, 4, 10),
         function(p) compound_symmetry(p)@param_params$link_rho@link_bounds,
         numeric(2)))
#>            [,1] [,2]
#> [1,] -1.0000000    1
#> [2,] -0.3333333    1
#> [3,] -0.1111111    1
```
