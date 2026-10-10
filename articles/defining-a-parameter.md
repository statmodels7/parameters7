# Defining a parameter

A covariance matrix is not a vector of free numbers. A map has to carry
an unconstrained vector, in which an optimizer can step, to a matrix
that is positive definite at every point of that vector, and it has to
carry the derivatives with it, because a fitting routine differentiates
through the map. A model may need a map that the package does not ship.
This vignette builds one, and shows at each step what the package
provides for it.

The example is an MA(1) autocovariance. For
$`x_k = e_k + \theta e_{k-1}`$ with $`\operatorname{Var}(e) = \sigma^2`$
the covariance is tridiagonal, $`\sigma^2(1 + \theta^2)`$ on the
diagonal and $`\sigma^2\theta`$ beside it, and it is positive definite
for every $`\lvert\theta\rvert < 1`$. The chart is therefore
$`(\log\sigma,\ \operatorname{atanh}\theta)`$, two unconstrained numbers
whatever the dimension.

## The minimum

A new parameter is a subclass of `parameter` (or of `matrix_parameter`,
if its value is a symmetric matrix and a consumer needs its determinant
and its solves) with a method for
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md).
The value carries the dimension labels `v1`, `v2`, … on both margins, as
the families of the package do.

``` r

Ma1Param <- S7::new_class("Ma1Param", parent = matrix_parameter)

ma1_band <- function(p, d, o) {
  m <- diag(d, p)
  if (p > 1) {
    i <- seq_len(p - 1)
    m[cbind(i, i + 1)] <- o
    m[cbind(i + 1, i)] <- o
  }
  dimnames(m) <- rep(list(paste0("v", seq_len(p))), 2)
  m
}

S7::method(param_value, Ma1Param) <- function(s, eta, ...) {
  p <- s@dimension
  s2 <- exp(2 * eta[[1]])
  th <- tanh(eta[[2]])
  ma1_band(p, s2 * (1 + th^2), s2 * th)
}
```

The constructor fills in what the class records about itself. `rank` and
`null_basis` are properties of the *family*, not of a point: they must
not move with the free vector, so a family whose null space would depend
on its coordinates is rejected at construction instead of recorded.

``` r

ma1 <- function(p, class = Ma1Param) {
  class(
    param_name = "ma1",
    dimension = as.integer(p),
    n_free = 2L,
    free_names = c("log_scale", "z_theta"),
    rank = as.integer(p),          # positive definite everywhere on the chart
    null_basis = matrix(0, p, 0),
    param_params = list()
  )
}

s <- ma1(5)
round(param_value(s, c(log(1.4), atanh(0.6))), 4)
#>        v1     v2     v3     v4     v5
#> v1 2.6656 1.1760 0.0000 0.0000 0.0000
#> v2 1.1760 2.6656 1.1760 0.0000 0.0000
#> v3 0.0000 1.1760 2.6656 1.1760 0.0000
#> v4 0.0000 0.0000 1.1760 2.6656 1.1760
#> v5 0.0000 0.0000 0.0000 1.1760 2.6656
```

Every generic except the inverse map
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
now has a method.

## The base-class methods

The derivatives, the log-determinant and its derivatives, the solve and
the factor all have methods on the base classes, so the family can be
used at once:

``` r

eta <- c(log(1.4), atanh(0.6))
str(param_d1(s, eta), max.level = 1)
#> List of 2
#>  $ log_scale: num [1:5, 1:5] 5.33 2.35 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ z_theta  : num [1:5, 1:5] 1.51 1.25 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
param_logdet(s, eta)
#> [1] 3.80883
round(param_solve(s, eta, b = c(1, 0, 0, 0, 0)), 6)
#>           [,1]
#> [1,]  0.508225
#> [2,] -0.301637
#> [3,]  0.175486
#> [4,] -0.096131
#> [5,]  0.042411
```

[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reports which quantities come from the base class: the derivatives of
the value are then finite differences, and the log-determinant and its
first two derivatives come from an eigendecomposition.

``` r

unlist(param_is_numerical(s))
#>       param_d1       param_d2       param_d3       param_d4   param_logdet 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#>  param_dlogdet param_d2logdet param_d3logdet param_d4logdet 
#>           TRUE           TRUE           TRUE           TRUE
```

For a family of the package every entry is `FALSE`, because every
quantity is written out:

``` r

unlist(param_is_numerical(log_cholesky(3)))
#>       param_d1       param_d2       param_d3       param_d4   param_logdet 
#>          FALSE          FALSE          FALSE          FALSE          FALSE 
#>  param_dlogdet param_d2logdet param_d3logdet param_d4logdet 
#>          FALSE          FALSE          FALSE          FALSE
```

The default second derivative and the default log-determinant
derivatives read the analytic first derivative when the family supplies
one, and the value otherwise; the default third and fourth derivatives
apply one stencil to the value. Supplying the first derivative therefore
improves the second derivative and the log-determinant derivatives.

## Adding a closed form

A registered method takes over through dispatch. The MA(1) derivatives
are elementary: with $`t = \tanh z`$ and $`t' = 1 - t^2`$,

``` math
\frac{\partial\Sigma}{\partial\log\sigma} = 2\Sigma,\qquad
  \frac{\partial\Sigma}{\partial z} = \sigma^2 t'
  \begin{pmatrix} 2t & 1 \\ 1 & 2t \end{pmatrix}\text{-banded.}
```

``` r

S7::method(param_d1, Ma1Param) <- function(s, eta, ...) {
  p <- s@dimension
  s2 <- exp(2 * eta[[1]])
  th <- tanh(eta[[2]])
  dt <- 1 - th^2
  stats::setNames(
    list(ma1_band(p, 2 * s2 * (1 + th^2), 2 * s2 * th),
         ma1_band(p, s2 * 2 * th * dt, s2 * dt)),
    s@free_names
  )
}

unlist(param_is_numerical(s))["param_d1"]
#> param_d1 
#>    FALSE
```

The [`setNames()`](https://rdrr.io/r/stats/setNames.html) is required:
the components of a derivative list are named by `free_names`, and those
of a second derivative by
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md).
A method that returns the right numbers in an unnamed list fails the
ninth check of the validator.

The default second derivative now differences this first derivative, and
the third and fourth orders remain stencils on the value. The following
chunk compares the first derivative with a central difference of the
value:

``` r

h <- 1e-6
fd <- lapply(1:2, function(j) {
  e1 <- eta; e1[j] <- e1[j] + h
  e0 <- eta; e0[j] <- e0[j] - h
  (param_value(s, e1) - param_value(s, e0)) / (2 * h)
})
max(abs(param_d1(s, eta)[[1]] - fd[[1]]))
#> [1] 3.192078e-10
max(abs(param_d1(s, eta)[[2]] - fd[[2]]))
#> [1] 2.024469e-10
```

## The log-determinant

A likelihood with a covariance contains $`\log\lvert\Sigma\rvert`$, so a
family can register methods for it and its derivatives. For a family
without such a method it is computed from an eigendecomposition; here
the matrix is a symmetric tridiagonal Toeplitz matrix, whose eigenvalues
are known in closed form,

``` math
\lambda_k = \sigma^2\bigl(1 + \theta^2 + 2\theta\cos\tfrac{k\pi}{p+1}\bigr),
  \qquad k = 1,\dots,p,
```

so the determinant costs $`p`$ cosines and no decomposition:

``` r

S7::method(param_logdet, Ma1Param) <- function(s, eta, ...) {
  p <- s@dimension
  s2 <- exp(2 * eta[[1]])
  th <- tanh(eta[[2]])
  sum(log(s2 * (1 + th^2 + 2 * th * cos(seq_len(p) * pi / (p + 1)))))
}

c(closed = param_logdet(s, eta),
  reference = as.numeric(determinant(param_value(s, eta),
                                     logarithm = TRUE)$modulus))
#>    closed reference 
#>   3.80883   3.80883
```

## Validating it

[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
checks that the value is symmetric and positive definite across the
chart, that the inverse map recovers the free vector, that the first and
second derivatives agree with a numerical reference, that the
log-determinant and its gradient and Hessian agree with the value, that
the solve and the factor reproduce the matrix, and that the shapes and
names are as declared.

``` r

res <- check_parameter(s, verbose = FALSE)
res
#>                check      status    statistic
#> 1         membership          OK 0.000000e+00
#> 2         round trip NOT CHECKED           NA
#> 3  first derivatives          OK 3.878954e-11
#> 4 second derivatives NOT CHECKED           NA
#> 5    log-determinant          OK 3.188277e-16
#> 6    logdet gradient NOT CHECKED           NA
#> 7     logdet hessian NOT CHECKED           NA
#> 8   solve and factor          OK 5.087223e-16
#> 9   shapes and names          OK           NA
```

A quantity that comes from the base class would be compared with a
numerical differentiation, which is the same arithmetic twice and would
agree however wrong the family is. Such a quantity is reported as **NOT
CHECKED**, so the report distinguishes what was verified from what was
only computed.

A validator is tested by breaking something and confirming that the
check fails:

``` r

Ma1Broken <- S7::new_class("Ma1Broken", parent = Ma1Param)

S7::method(param_d1, Ma1Broken) <- function(s, eta, ...) {
  d <- S7::method(param_d1, Ma1Param)(s, eta)   # correct, then spoiled
  d[[1]] <- d[[1]] * 1.05                       # 5 percent too large
  d
}

bad <- check_parameter(ma1(5, class = Ma1Broken), verbose = FALSE)
bad[bad$status != "OK", ]
#>                check      status statistic
#> 2         round trip NOT CHECKED        NA
#> 3  first derivatives        FAIL      0.05
#> 4 second derivatives NOT CHECKED        NA
#> 6    logdet gradient NOT CHECKED        NA
#> 7     logdet hessian NOT CHECKED        NA
```

An error of one part in twenty is caught, and the checks that were not
checked before remain NOT CHECKED.

## Rank, the null space and the solve

A family may be rank deficient, and then some results change.
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
on a singular matrix records the rank, reports the **log
pseudo-determinant**, the determinant itself being $`-\infty`$, and
rejects the solve:

``` r

sm <- scaled_matrix(matrix(1, 2, 2))       # rank 1 by construction
c(dimension = sm@dimension, rank = sm@rank)
#> dimension      rank 
#>         2         1
param_logdet(sm, 0)                        # log of the one non-zero eigenvalue
#> [1] 0.6931472
try(param_solve(sm, 0))
#> Error : 'scaled' is rank deficient (1 of 2), so it has no inverse. A consumer of
#>   an improper prior needs the quadratic form and the log
#>   pseudo-determinant, not a pseudo-inverse: assemble the matrix the
#>   model actually inverts and solve that.
```

The rejection is deliberate. A consumer of an improper prior needs the
quadratic form and the log pseudo-determinant, and the matrix that it
inverts is $`X'X + \lambda P`$, which is non-singular. A pseudo-inverse
returned by
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
would be a different quantity from the inverse.

## Names of coordinates and quantities

A free value is unconstrained, so a name that suggests a bounded
quantity for a free value would mislead. `free_names` records the
*transform*, and
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
reports the quantities that the family describes, with the Jacobian from
the free vector for the delta method:

``` r

ar <- autoregressive(4, order = 2)
ar@free_names
#> [1] "log_scale" "z_pacf1"   "z_pacf2"
rd <- param_readable(ar, c(0, 0.5, 0.2))
data.frame(quantity = rownames(rd$jacobian), value = signif(rd$value, 6),
           row.names = NULL)
#>   quantity    value
#> 1    scale 1.000000
#> 2    pacf1 0.462117
#> 3    pacf2 0.197375
#> 4     phi1 0.370907
#> 5     phi2 0.197375
```

The coordinate `z_pacf1` is the inverse hyperbolic tangent of the first
partial autocorrelation, `pacf1`, and `phi1` is the first autoregressive
coefficient, which does not appear among the coordinates. The intervals
of the coefficients are built on the identity scale, the stationary
region not being a box.

## Composing

Five wrappers build a parameter from others:
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md),
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md),
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
and
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md).
Each assembles its derivatives from those of its parts, with no new
derivation:

``` r

blk <- block_diag(list(ma1(3), log_cholesky(2)))
c(n_free = blk@n_free, dimension = blk@dimension)
#>    n_free dimension 
#>         5         5
blk@free_names
#> [1] "b1_log_scale" "b1_z_theta"   "b2_log_L1"    "b2_log_L2"    "b2_L2.1"
```

In a block diagonal every component whose indices span two blocks is
exactly zero at every order. `kron_identity(s, m)` repeats one block
$`m`$ times over one shared free vector, as grouped random effects need;
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
writes a covariance as $`DRD`$, so the coordinates are the linked
standard deviations followed by those of the correlation block;
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
is a non-negative combination of fixed matrices, and its log-determinant
derivatives are the cyclic trace expansion;
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
gives the family whose value is the inverse of another family’s.

[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
reads its rank from the components **stacked and individually
normalized**, never from the assembled matrix. The null space of a sum
of positive semidefinite matrices is the intersection of theirs and does
not move with the weights, while a count taken from the assembled matrix
falls as the weights spread apart.

## Summary

- **Minimum to define a parameter:** a subclass of `parameter` or
  `matrix_parameter` and a
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  method, with the constructor recording `n_free`, `free_names`, and for
  a matrix the `dimension`, `rank` and `null_basis`.
- Rank and null space are properties of the family and must not move
  with the free vector; a family whose null space would move is
  rejected.
- Every generic except
  [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
  has a method on the base classes, and
  [`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
  reports which quantities use it.
- Closed forms can be registered one at a time, each taking over through
  dispatch; a closed first derivative improves the default second
  derivative and the log-determinant derivatives.
- [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  verifies what can be verified and reports the rest as NOT CHECKED.
- `free_names` names the transform;
  [`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
  names the quantity and carries the Jacobian.
- The five composition wrappers assemble their derivatives from those of
  their parts.
