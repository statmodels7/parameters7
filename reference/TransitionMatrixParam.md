# Transition Matrix Parameter

The S7 class of row-stochastic matrices, each row on the open simplex in
the additive log-ratio parametrization. Its value is a matrix but **not
a symmetric** one, so it inherits
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
directly, never
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md):
there is no `rank`, no `null_basis` and no `role`, and a transition
matrix has no log-determinant, solve or factor to be asked for.

[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
builds one. The rows are independent in the parametrization, so every
derivative array is block diagonal by row and the family reuses
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)'s
kernels row by row.

## Usage

``` r
TransitionMatrixParam(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list()
)
```

## Arguments

- param_name:

  A single character string naming the family, used in the error
  messages the validators raise and in the object's `print` output.
  `"log_cholesky"`, `"ar1"` and so on.

- n_free:

  The length \\d\\ of the free vector: a single non-negative integer.
  Zero is legal and describes a family with nothing to estimate. The
  validator rejects a vector, a negative value, or a length that
  disagrees with `free_names`.

- free_names:

  A character vector of length `n_free`, one label per free value, in
  the order the free vector holds them. Fixed at construction and part
  of the interface: consumers build their parameter tables from these
  labels, so the ordering is not free to change. The validator rejects a
  duplicated label and a length other than `n_free`.

- param_params:

  A list of whatever the family needs in order to evaluate itself, read
  only by that family's own methods.
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  stores the row and column index of each free value here;
  [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
  stores its component matrices. Nothing outside the family looks inside
  it.

## Value

An object of class `TransitionMatrixParam`, a subclass of
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
adding no properties of its own: `param_name` is `"transition_matrix"`,
`n_free` is \\K(K-1)\\, `free_names` is `alr{i}.{j}` row by row, and
`param_params` holds `n_state`.

## See also

[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
the constructor,
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md),
which is one row of this, and
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
for the properties this inherits and the generics it does not get.

## Examples

``` r
# A matrix, but not a symmetric one, so not a matrix_parameter.
s <- transition_matrix(3)
c(parameter = S7::S7_inherits(s, parameter),
  matrix_parameter = S7::S7_inherits(s, matrix_parameter))
#>        parameter matrix_parameter 
#>             TRUE            FALSE 

# K(K-1) free values: one simplex per row.
rbind(K = 2:5, n_free = vapply(2:5, function(k)
  transition_matrix(k)@n_free, integer(1)))
#>        [,1] [,2] [,3] [,4]
#> K         2    3    4    5
#> n_free    2    6   12   20
```
