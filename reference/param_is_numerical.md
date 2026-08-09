# Which of a Parameter's Quantities Come From the Base Class

Reports, for each of the five derivative quantities, whether the
parameter supplies its own method or falls back to the one registered on
[`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

## Usage

``` r
param_is_numerical(s)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

## Value

A named logical vector over `param_d1`, `param_d2`, `param_logdet`,
`param_dlogdet` and `param_d2logdet`, `TRUE` where the base-class method
is in force.

## Details

The distinction that matters is whether an independent check exists. A
derivative computed by finite differences cannot be checked against a
finite difference, and a log-determinant read off an eigendecomposition
cannot be checked against an eigendecomposition; the comparison is the
same arithmetic twice, and it agrees however wrong the parameter is.
[`check_parameter`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
uses this to report such a quantity as not checked rather than as
passed.

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
and
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
are deliberately absent. Their base-class versions are a Cholesky
factorization, which is exact whoever performs it, and the validator
compares them with [`base::solve`](https://rdrr.io/r/base/solve.html)
either way. Calling them numerical would suggest an approximation that
is not there.

## See also

[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md),
[`param_free`](https://statmodels7.github.io/parameters7/reference/param_free.md)

## Examples

``` r
# everything closed form except the second-order log-determinant of a
# parameter that does not supply one
param_is_numerical(log_cholesky(3))
#>       param_d1       param_d2       param_d3       param_d4   param_logdet 
#>          FALSE          FALSE          FALSE          FALSE          FALSE 
#>  param_dlogdet param_d2logdet param_d3logdet param_d4logdet 
#>          FALSE          FALSE          FALSE          FALSE 
param_is_numerical(diagonal_matrix(3))
#>       param_d1       param_d2       param_d3       param_d4   param_logdet 
#>          FALSE          FALSE          FALSE          FALSE          FALSE 
#>  param_dlogdet param_d2logdet param_d3logdet param_d4logdet 
#>          FALSE          FALSE          FALSE          FALSE 
```
