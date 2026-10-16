# Simplex Parameter

`SimplexParam` is the S7 class of probability vectors on the open
simplex, in the additive log-ratio parametrization. It inherits
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
directly and not
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
because its value is a vector. The class therefore lacks the properties
`dimension`, `rank` and `null_basis`, and the generics
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
are not defined for it, so a call for the log-determinant of a
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
  messages that the validators raise and in the object's `print` output.
  Examples are `"log_cholesky"` and `"ar1"`.

- n_free:

  The length \\d\\ of the free vector: a single non-negative integer.
  Zero is allowed and describes a family with nothing to estimate. The
  validator rejects a vector, a negative value, or a length that
  disagrees with `free_names`.

- free_names:

  A character vector of length `n_free`, one label per free value, in
  the order in which the free vector holds them. Fixed at construction
  and part of the interface: consumers build their parameter tables from
  these labels, so the ordering is not free to change. The validator
  rejects a duplicated label and a length other than `n_free`.

- param_params:

  A list of whatever the family needs in order to evaluate itself, read
  only by that family's own methods.
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  stores the row and column index of each free value here;
  [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
  stores its component matrices. Only the methods of that family read
  the list.

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
which is one of these per row,
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
for the properties that this class inherits, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the generics that are not defined for it.

## Examples

``` r
# A probability vector is not a matrix, so it inherits parameter() alone.
s <- simplex(3)
c(parameter = S7::S7_inherits(s, parameter),
  matrix_parameter = S7::S7_inherits(s, matrix_parameter))
#>        parameter matrix_parameter 
#>             TRUE            FALSE 

# A simplex of K categories has K - 1 free values.
vapply(2:5, function(k) simplex(k)@n_free, integer(1))
#> [1] 1 2 3 4

# The log-determinant is not defined for a simplex.
try(param_logdet(s, c(0, 0)))
#> Error : Can't find method for `param_logdet(<parameters7::SimplexParam>)`.
```
