# Validate a Covariance Parameter

Runs a battery of numerical checks on a parametrization, each against a
route the implementation does not itself take, and reports what passed,
what failed and what could not be checked. Write a family of your own,
call this on it, and read the table: a `FAIL` names the quantity whose
closed form is wrong, and the statistic beside it is the relative size
of the disagreement.

It is also the validator the package holds itself to. Measured over all
fifteen constructors: thirteen matrix families pass all nine checks and
the two that are not matrices pass all seven, with nothing skipped
anywhere.

## Usage

``` r
check_parameter(s, tol = 1e-06, verbose = TRUE)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md).
  Anything else throws `'s' must inherit from class 'parameter'.`

- tol:

  The relative tolerance a check has to meet, defaulting to `1e-6`. The
  comparisons are relative to the larger of 1 and the size of the
  reference, so the tolerance is dimensionless. `1e-6` is loose enough
  for the third and fourth derivative comparisons, which rest on
  stencils good to about \\10^{-5}\\, and tight enough to catch an error
  of one part in a thousand.

- verbose:

  Whether to print the table, `TRUE` by default. The value is returned
  either way.

## Value

Invisibly, a data frame with one row per check. For a
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
it has nine rows and the columns

- `check`:

  character, the names listed above.

- `status`:

  character, `"OK"`, `"FAIL"` or `"NOT CHECKED"`.

- `statistic`:

  numeric, the worst relative discrepancy, `NA` for a check that was
  skipped or has no number.

For a family that is not a matrix it has seven rows and the columns
`check`, `status` and `note`, the last being **character** and holding a
formatted number or the empty string. Test the `status` column, which is
common to both.

## Not checked is a third verdict

A quantity that comes from a numerical fallback is reported as **NOT
CHECKED**, never as passed. Comparing a finite difference against a
finite difference is the same arithmetic twice, and it agrees however
wrong the parametrization is.
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
is what decides, and a family supplying only
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
therefore comes back with six of the nine rows skipped: it has nothing
an independent route could contradict.

## The nine checks, and the route each is held against

1.  **membership**: the matrix is symmetric and positive semidefinite, a
    full-rank family has a strictly positive smallest eigenvalue, and a
    rank-deficient one annihilates its declared null space. The last is
    tested through the null basis, never by counting eigenvalues, a
    count of small eigenvalues not being scale invariant.

2.  **round trip**:
    [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
    recovers the free vector from the matrix, where the family
    implements an inverse.

3.  **first derivatives** against one central difference of
    [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md).

4.  **second derivatives** against one central difference of the
    analytic first derivatives.

5.  **log-determinant** against the sum of the logs of the eigenvalues
    the declared rank keeps.

6.  **logdet gradient** against \\\mathrm{tr}(M^{+} \partial_k M)\\,
    with the pseudo-inverse formed from an eigendecomposition instead of
    from
    [`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md),
    so the two routes share no arithmetic.

7.  **logdet hessian** against one central difference of the analytic
    gradient.

8.  **solve and factor** against
    [`base::solve()`](https://rdrr.io/r/base/solve.html), where the
    family is of full rank; skipped for a deficient one, which has no
    inverse.

9.  **shapes and names**: the declared dimension, the lengths and the
    names match what the methods return. Structural, so it carries no
    statistic.

Each check reports the worst discrepancy over the free vectors
[`sweep_etas()`](https://statmodels7.github.io/parameters7/reference/sweep_etas.md)
supplies, four of which are drawn at random.

## Neither branch touches the caller's random stream

Both batteries draw from a fixed seed, so two calls on the same family
report the same statistics, and both put back the `.Random.seed` they
found on entry, so a call in the middle of a simulation leaves that
simulation unchanged. The matrix branch seeds once before
[`sweep_etas()`](https://statmodels7.github.io/parameters7/reference/sweep_etas.md);
the branch for a family that is not a matrix seeds each of its three
draws. Restoring is
[`capture_seed()`](https://statmodels7.github.io/parameters7/reference/capture_seed.md)
and
[`restore_seed()`](https://statmodels7.github.io/parameters7/reference/capture_seed.md),
through [`base::on.exit()`](https://rdrr.io/r/base/on.exit.html), so it
happens even when a check signals.

## A family that is not a matrix gets a different battery

[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
have no log-determinant, no solve and no factor, so seven other checks
run instead: the inverse round trip, the four derivative orders against
the single-stencil construction, that the value stays on the simplex,
and that every derivative component sums to zero over the value index.
**The returned table has different columns in that case**; see
**Value**.

## See also

[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md),
which decides what is checkable, and
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
for the null space check 1 uses.

## Examples

``` r
set.seed(1)

# A family that passes everything.
r <- check_parameter(log_cholesky(3))
#> Parameter: log_cholesky   (3 x 3, rank 3, 6 free)
#>   [OK         ] membership           0.00e+00
#>   [OK         ] round trip           3.25e-16
#>   [OK         ] first derivatives    2.90e-11
#>   [OK         ] second derivatives   3.67e-11
#>   [OK         ] log-determinant      6.13e-15
#>   [OK         ] logdet gradient      1.61e-13
#>   [OK         ] logdet hessian       0.00e+00
#>   [OK         ] solve and factor     6.04e-15
#>   [OK         ] shapes and names  
#>   9 passed, 0 failed, 0 not checked
r$status
#> [1] "OK" "OK" "OK" "OK" "OK" "OK" "OK" "OK" "OK"

# A rank-deficient penalty passes the same battery, with the solve skipped:
# it has no inverse, and the null-space check replaces the definiteness one.
P <- crossprod(diff(diag(6), differences = 2))
d <- check_parameter(scaled_matrix(P))
#> Parameter: scaled   (6 x 6, rank 4, 1 free)
#>   [OK         ] membership           3.19e-16
#>   [OK         ] round trip           5.55e-17
#>   [OK         ] first derivatives    2.40e-11
#>   [OK         ] second derivatives   2.40e-11
#>   [OK         ] log-determinant      4.33e-15
#>   [OK         ] logdet gradient      1.78e-15
#>   [OK         ] logdet hessian       0.00e+00
#>   [NOT CHECKED] solve             
#>   [OK         ] shapes and names  
#>   8 passed, 0 failed, 1 not checked
d[d$status != "OK", ]
#>   check      status statistic
#> 8 solve NOT CHECKED        NA

# A family that is not a matrix gets the seven-row battery and a `note`
# column in place of `statistic`.
v <- check_parameter(simplex(4), verbose = FALSE)
names(v)
#> [1] "check"  "status" "note"  
nrow(v)
#> [1] 7

# What a real defect looks like. This first derivative is 5 per cent wrong.
Wrong <- S7::new_class("Wrong", parent = matrix_parameter)
S7::method(param_value, Wrong) <- function(s, eta, ...) {
  m <- diag(rep(exp(eta[1]), 2))
  dimnames(m) <- list(c("v1", "v2"), c("v1", "v2"))
  m
}
S7::method(param_d1, Wrong) <- function(s, eta, ...) {
  list(log_s = 1.05 * diag(rep(exp(eta[1]), 2)))
}
w <- Wrong(param_name = "wrong", n_free = 1L, free_names = "log_s",
           param_params = list(), dimension = 2L, rank = 2L,
           null_basis = matrix(numeric(0), 2, 0), role = "either")
bad <- check_parameter(w, verbose = FALSE)
bad[bad$status == "FAIL", ]
#>               check status statistic
#> 3 first derivatives   FAIL      0.05
#> 9  shapes and names   FAIL        NA
```
