# Simplex Parameter

The S7 class of probability vectors on the open simplex, in the additive
log-ratio parametrization. It inherits
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
**directly**, never
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
its value being a vector: there is no `dimension`, no `rank`, no
`null_basis` and no `role`, and
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
have no method for it, so asking for the log-determinant of a
probability vector fails at dispatch.

[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
builds one. The four derivative orders are all closed form.

## Usage

``` r
SimplexParam(
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

An object of class `SimplexParam`, a subclass of
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
adding no properties of its own: `param_name` is `"simplex"`, `n_free`
is \\K - 1\\, `free_names` is `alr1` ... `alr(K-1)`, and `param_params`
holds `n_cat`.

## See also

[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md),
the constructor,
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
which is one of these per row, and
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
for the properties this inherits and the four generics it does not get.

## Examples

``` r
# A probability vector is not a matrix, so it inherits parameter() alone.
s <- simplex(3)
c(parameter = S7::S7_inherits(s, parameter),
  matrix_parameter = S7::S7_inherits(s, matrix_parameter))
#>        parameter matrix_parameter 
#>             TRUE            FALSE 

# K - 1 free values for K categories.
vapply(2:5, function(k) simplex(k)@n_free, integer(1))
#> [1] 1 2 3 4

# And no log-determinant to ask for.
try(param_logdet(s, c(0, 0)))
#> Error : Can't find method for `param_logdet(<parameters7::SimplexParam>)`.
```
