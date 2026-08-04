# Validate a Covariance Parameter

Runs a battery of numerical checks on a parameter, each against a route
the implementation does not itself take, and reports what passed, what
failed and what could not be checked.

## Usage

``` r
check_parameter(s, tol = 1e-06, verbose = TRUE)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- tol:

  The relative tolerance for the comparisons.

- verbose:

  Whether to print the table.

## Value

Invisibly, a data frame with columns `check`, `status` and `statistic`.

## Details

A quantity that comes from a numerical fallback is reported as **not
checked** rather than as passed. Comparing a finite difference against a
finite difference is the same arithmetic twice, and it agrees however
wrong the parameter is; saying so is the difference between a validator
and a formality.

The checks are:

1.  **membership**: the matrix is symmetric and positive semidefinite, a
    full-rank family has a positive smallest eigenvalue, and a
    rank-deficient one annihilates its declared null space. The last is
    tested through the null basis rather than by counting eigenvalues,
    because a count is not scale invariant.

2.  **round trip**:
    [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
    recovers the free vector from the matrix, where the family
    implements it.

3.  **first derivatives** against one central difference of
    [`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md).

4.  **second derivatives** against one central difference of the
    analytic first derivatives.

5.  **log-determinant** against the sum of the logs of the eigenvalues
    the rank keeps.

6.  **log-determinant gradient** against \\\mathrm{tr}(M^{+} \partial_k
    M)\\, with the pseudo-inverse formed from an eigendecomposition
    rather than from
    [`param_solve`](https://statmodels7.github.io/parameters7/reference/param_solve.md).

7.  **log-determinant Hessian** against one central difference of the
    analytic gradient.

8.  **solve** against
    [`base::solve`](https://rdrr.io/r/base/solve.html), where the family
    is of full rank.

9.  **shapes**: the declared dimension, length and names match what the
    methods return.

## See also

[`param_is_numerical`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)

## Examples

``` r
invisible(check_parameter(log_cholesky(3)))
#> Parameter: log_cholesky   (3 x 3, rank 3, 6 free)
#>   [OK         ] membership           0.00e+00
#>   [OK         ] round trip           4.94e-16
#>   [OK         ] first derivatives    2.90e-11
#>   [OK         ] second derivatives   3.67e-11
#>   [OK         ] log-determinant      6.38e-15
#>   [OK         ] logdet gradient      2.30e-13
#>   [OK         ] logdet hessian       0.00e+00
#>   [OK         ] solve and factor     3.76e-15
#>   [OK         ] shapes and names  
#>   9 passed, 0 failed, 0 not checked

# a rank-deficient penalty passes the same battery
invisible(check_parameter(scaled_matrix(crossprod(diff(diag(6), differences = 2)))))
#> Parameter: scaled   (6 x 6, rank 4, 1 free)
#>   [OK         ] membership           3.29e-16
#>   [OK         ] round trip           1.11e-16
#>   [OK         ] first derivatives    2.40e-11
#>   [OK         ] second derivatives   2.40e-11
#>   [OK         ] log-determinant      1.33e-15
#>   [OK         ] logdet gradient      1.55e-15
#>   [OK         ] logdet hessian       0.00e+00
#>   [NOT CHECKED] solve             
#>   [OK         ] shapes and names  
#>   8 passed, 0 failed, 1 not checked
```
