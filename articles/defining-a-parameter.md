# Defining a parameter

A covariance matrix is not a vector of free numbers. Something has to
carry it from an unconstrained vector an optimizer can step in to a
matrix that is positive definite at every point of that vector, and it
has to carry the derivatives with it, because a fitting routine
differentiates through the map. A family that only worked for the
fourteen charts this package ships would have solved nothing. This
vignette builds one it does not ship, and shows what the package
guarantees about it at each step.

The example is an MA(1) autocovariance. For
$`x_k = e_k + \theta e_{k-1}`$ with $`\operatorname{Var}(e) = \sigma^2`$
the covariance is tridiagonal, $`\sigma^2(1 + \theta^2)`$ on the
diagonal and $`\sigma^2\theta`$ beside it, and it is positive definite
for every $`\lvert\theta\rvert < 1`$. The chart is therefore
$`(\log\sigma,\ \operatorname{atanh}\theta)`$, two unconstrained numbers
whatever the dimension.

## The minimum

A new parameter is a subclass of `parameter` (or of `matrix_parameter`,
if it is a matrix and a consumer will want its determinant and its
solves) and a method for
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md).
Nothing else is required.

``` r

Ma1Param <- S7::new_class("Ma1Param", parent = matrix_parameter)

ma1_band <- function(p, d, o) {
  m <- diag(d, p)
  if (p > 1) {
    i <- seq_len(p - 1)
    m[cbind(i, i + 1)] <- o
    m[cbind(i + 1, i)] <- o
  }
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
not move with the free vector, which is why a family whose null space
depends on its coordinates has to reject rather than record one.

``` r

ma1 <- function(p, class = Ma1Param) {
  class(
    param_name = "ma1",
    dimension = as.integer(p),
    n_free = 2L,
    free_names = c("log_scale", "z_theta"),
    rank = as.integer(p),          # positive definite everywhere on the chart
    null_basis = matrix(0, p, 0),
    role = "covariance",
    param_params = list()
  )
}

s <- ma1(5)
round(param_value(s, c(log(1.4), atanh(0.6))), 4)
#>        [,1]   [,2]   [,3]   [,4]   [,5]
#> [1,] 2.6656 1.1760 0.0000 0.0000 0.0000
#> [2,] 1.1760 2.6656 1.1760 0.0000 0.0000
#> [3,] 0.0000 1.1760 2.6656 1.1760 0.0000
#> [4,] 0.0000 0.0000 1.1760 2.6656 1.1760
#> [5,] 0.0000 0.0000 0.0000 1.1760 2.6656
```

That is the whole of the compulsory part. Every other generic now
answers.

## What comes free, and what it costs

The derivatives, the log-determinant and its derivatives, the solve and
the factor all have numerical fallbacks on the base class, so the family
is usable at once:

``` r

eta <- c(log(1.4), atanh(0.6))
str(param_d1(s, eta), max.level = 1)
#> List of 2
#>  $ log_scale: num [1:5, 1:5] 5.33 2.35 0 0 0 ...
#>  $ z_theta  : num [1:5, 1:5] 1.51 1.25 0 0 0 ...
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
says which of them are differences and which are formulas, so a reader
can tell what the object is doing:

``` r

unlist(param_is_numerical(s))
#>       param_d1       param_d2       param_d3       param_d4   param_logdet 
#>           TRUE           TRUE           TRUE           TRUE           TRUE 
#>  param_dlogdet param_d2logdet param_d3logdet param_d4logdet 
#>           TRUE           TRUE           TRUE           TRUE
```

Compare a shipped family, where the answer is empty because every order
is written out:

``` r

unlist(param_is_numerical(log_cholesky(3)))
#>       param_d1       param_d2       param_d3       param_d4   param_logdet 
#>          FALSE          FALSE          FALSE          FALSE          FALSE 
#>  param_dlogdet param_d2logdet param_d3logdet param_d4logdet 
#>          FALSE          FALSE          FALSE          FALSE
```

The cost is accuracy and time. A fallback of order $`k`$ applies **one**
stencil to the highest analytic order the family supplies, never a chain
of first differences, so supplying the first derivative improves every
order above it.

## Adding a closed form

Register one method and it takes over through dispatch. The MA(1)
derivatives are elementary: with $`t = \tanh z`$ and $`t' = 1 - t^2`$,

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

unlist(param_is_numerical(s))["d1"]
#> <NA> 
#>   NA
```

The [`setNames()`](https://rdrr.io/r/stats/setNames.html) is part of the
contract, no decoration: the components of a derivative list are named
by `free_names`, and those of a second derivative by
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md).
A method returning the right numbers in an unnamed list is caught by the
validator’s ninth check, which is how this one was written correctly on
the second attempt.

The higher orders are still numerical, and are now one stencil on this,
no longer on the value:

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

The contract asks for $`\log\lvert\Sigma\rvert`$ and its derivatives,
because a likelihood carrying a covariance carries that term. For most
families it is computed from a factorization; here the matrix is a
symmetric tridiagonal Toeplitz, whose eigenvalues are known in closed
form,

``` math
\lambda_k = \sigma^2\bigl(1 + \theta^2 + 2\theta\cos\tfrac{k\pi}{p+1}\bigr),
  \qquad k = 1,\dots,p,
```

so the determinant costs $`p`$ cosines and no decomposition at all:

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
asks what can be asked: that the value is symmetric and positive
definite across the chart, that every derivative order agrees with a
numerical reference, that the log-determinant and its derivatives agree
with the value’s, and that the solve and the factor reproduce the
matrix.

``` r

res <- check_parameter(s, verbose = FALSE)
res
#>                check      status    statistic
#> 1         membership          OK 0.000000e+00
#> 2         round trip NOT CHECKED           NA
#> 3  first derivatives          OK 2.813683e-11
#> 4 second derivatives NOT CHECKED           NA
#> 5    log-determinant          OK 3.188277e-16
#> 6    logdet gradient NOT CHECKED           NA
#> 7     logdet hessian NOT CHECKED           NA
#> 8   solve and factor          OK 5.087223e-16
#> 9   shapes and names        FAIL           NA
```

An order that comes from a fallback is compared against a numerical
differentiation of the order below it, which is the same arithmetic
twice and would agree however wrong the family is. Those are reported as
**unchecked**, never as passed, so the report distinguishes what was
verified from what was merely computed.

The way to test a validator is to break something and confirm it
complains:

``` r

Ma1Broken <- S7::new_class("Ma1Broken", parent = Ma1Param)

S7::method(param_d1, Ma1Broken) <- function(s, eta, ...) {
  d <- S7::method(param_d1, Ma1Param)(s, eta)   # correct, then spoiled
  d[[1]] <- d[[1]] * 1.05                       # five per cent out
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
#> 9   shapes and names        FAIL        NA
```

An error of one part in twenty is caught, and the checks that were
already unchecked stay so.

## Rank, the null space and the solve

A family is allowed to be rank deficient, and then the honest answers
change.
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
on a singular matrix records the rank, reports the **log
pseudo-determinant**, the determinant itself being $`-\infty`$, and
refuses the solve:

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

The refusal is deliberate. What a consumer of an improper prior needs is
the quadratic form and the log pseudo-determinant; the matrix it
actually inverts is $`X'X + \lambda P`$, which is non-singular. A
pseudo-inverse returned under the name of a solve would be a different
object with the same name.

## A name says the coordinate, not the quantity

A free value is unconstrained by construction, so a name that promises a
bounded quantity and reports a free one misleads. `free_names` records
the *transform*, and
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
reports what the family is about, with the Jacobian from the free vector
so a delta method can follow:

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

The two are different numbers on purpose: the coordinate is a partial
autocorrelation and `phi1` is the first autoregressive coefficient,
which appears nowhere in $`\Sigma`$’s coordinates. The coefficients are
intervalled on the identity scale, the stationary region not being a
box.

## Composing

Four wrappers build a parameter from others, and none of them rederives
anything, because the free values of one part do not enter another:

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
$`m`$ times over one free vector, as grouped random effects need;
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
writes a covariance as $`DRD`$, so the coordinates are standard
deviations and correlations;
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
is a non-negative combination of fixed matrices, and its log-determinant
derivatives are the cyclic trace expansion.

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
  with the free vector; a family whose null space does move rejects
  instead.
- Everything else has a numerical fallback, and
  [`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
  says which is which.
- Register closed forms one at a time; each takes over through dispatch,
  and each improves the fallbacks above it.
- [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  verifies what can be verified and reports the rest as unchecked.
- `free_names` names the transform;
  [`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
  names the quantity and carries the Jacobian.
- The four composition wrappers rederive nothing.
