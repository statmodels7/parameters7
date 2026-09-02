# Which of a Parameter's Quantities Come From the Base Class

Reports, one derivative quantity at a time, whether the parameter has a
method of its own or takes the numerical one registered on
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md).
Ask it before trusting a fourth-order derivative of a family you did not
write: `FALSE` everywhere means every quantity is a closed form, and
`TRUE` somewhere marks a component whose accuracy is a stencil's.

## Usage

``` r
param_is_numerical(s)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md).

## Value

A named logical vector, `TRUE` where the base-class method is in force.
Its length depends on the branch: **nine** entries for a
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
over `param_d1`, `param_d2`, `param_d3`, `param_d4`, `param_logdet`,
`param_dlogdet`, `param_d2logdet`, `param_d3logdet` and
`param_d4logdet`; **four** for a family that is not a matrix, over the
derivative orders alone, the log-determinant not existing there.

## Why the question is worth asking

Whether an independent check exists depends on the answer. A derivative
computed by finite differences cannot be checked against a finite
difference, and a log-determinant read off an eigendecomposition cannot
be checked against an eigendecomposition: the comparison is the same
arithmetic twice, and it agrees however wrong the parametrization is.
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
calls this and reports such a quantity as **not checked** instead of as
passed.

## What is not in the list

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
are absent deliberately. Their base-class versions are a Cholesky
factorization, which is exact whoever performs it, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
compares them against
[`base::solve()`](https://rdrr.io/r/base/solve.html) either way. Listing
them as numerical would suggest an approximation that is not there.
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
is absent because a family that does not implement it has nothing at
all, and
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
because its base method throws.

## Every shipped family answers FALSE

Measured over all fifteen constructors in this package, at every
component: none of them uses a numerical route. The fallbacks exist for
a family written elsewhere, and the example below builds one to show
what a `TRUE` looks like.

## See also

[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
which reports a numerical component as not checked, and
[`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md)
and its higher-order siblings, the routes a `TRUE` names.

## Examples

``` r
# Every family in the package is closed form throughout.
param_is_numerical(log_cholesky(3))
#>       param_d1       param_d2       param_d3       param_d4   param_logdet 
#>          FALSE          FALSE          FALSE          FALSE          FALSE 
#>  param_dlogdet param_d2logdet param_d3logdet param_d4logdet 
#>          FALSE          FALSE          FALSE          FALSE 
any(param_is_numerical(diagonal_matrix(3)))
#> [1] FALSE

# Nine components for a matrix family, four for one that is not.
c(matrix = length(param_is_numerical(log_cholesky(3))),
  not_matrix = length(param_is_numerical(simplex(3))))
#>     matrix not_matrix 
#>          9          4 

# A family written from scratch, with param_value() alone: every derivative
# order is then numerical, and so is the whole log-determinant block.
Toy <- S7::new_class("Toy", parent = matrix_parameter)
S7::method(param_value, Toy) <- function(s, eta, ...) {
  m <- diag(rep(exp(eta[1]), 2))
  dimnames(m) <- list(c("v1", "v2"), c("v1", "v2"))
  m
}
toy <- Toy(param_name = "toy", n_free = 1L, free_names = "log_s",
           param_params = list(), dimension = 2L, rank = 2L,
           null_basis = matrix(numeric(0), 2, 0))
param_is_numerical(toy)
#>       param_d1       param_d2       param_d3       param_d4   param_logdet 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#>  param_dlogdet param_d2logdet param_d3logdet param_d4logdet 
#>           TRUE           TRUE           TRUE           TRUE 

# The numerical route is usable, and here it can be checked by hand: the
# matrix is exp(eta) times the identity, so every derivative is itself.
c(numerical = param_d1(toy, 0.4)[[1]][1, 1], exact = exp(0.4))
#> numerical     exact 
#>  1.491825  1.491825 
```
