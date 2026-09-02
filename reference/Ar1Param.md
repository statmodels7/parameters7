# AR(1) Parameter

The S7 class of first-order autoregressive covariance matrices: equal
variances and a correlation falling geometrically with the lag,
\\M\_{ij} = \sigma^2 \rho^{\|i-j\|}\\. Two free values at every
dimension, and the correlation is unrestricted in \\(-1, 1)\\, unlike
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)'s.

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
builds one. It shares its whole derivative skeleton with
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
through
[`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md),
[`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md)
and
[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md);
what differs is the pattern and the log-determinant terms.

## Usage

``` r
Ar1Param(
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

An object of class `Ar1Param`, a subclass of
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
adding no properties of its own. `param_params` holds `link_scale` and
`link_rho`, the second a
[`linkfunctions7::rhobit_link()`](https://statmodels7.github.io/linkfunctions7/reference/rhobit_link.html)
at every dimension. `n_free` is 2 and `rank` is \\p\\.

## See also

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
the constructor,
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
for the sibling family,
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for higher orders, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the properties this inherits.

## Examples

``` r
# Two free values whatever the dimension, and a rhobit correlation.
s <- ar1(4)
S7::S7_inherits(s, Ar1Param)
#> [1] TRUE
c(p4 = ar1(4)@n_free, p20 = ar1(20)@n_free)
#>  p4 p20 
#>   2   2 
s@free_names
#> [1] "log_scale" "z_rho"    

# Unlike compound symmetry, the bound does not move with the dimension.
t(vapply(c(2, 4, 20), function(p) ar1(p)@param_params$link_rho@link_bounds,
         numeric(2)))
#>      [,1] [,2]
#> [1,]   -1    1
#> [2,]   -1    1
#> [3,]   -1    1
```
