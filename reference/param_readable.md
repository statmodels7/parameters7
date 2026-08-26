# Quantities a Family Is About

Asks a family which interpretable quantities it stands for, and returns
them with everything a consumer needs to report them: the values, the
Jacobian of the map from the free vector, the scale each interval should
be built on, and a label for the block. It exists because **the free
vector is what a fit estimates, never what a reader reads**: nobody
reads the hyperbolic arc tangent of a partial autocorrelation, and the
quantity behind it cannot be recovered from a printed covariance either.

Use it to turn an estimate and its variance matrix into a table a reader
can use. The base class declares nothing, so a family whose matrix is
all there is to say returns `NULL` and a consumer falls back on
reporting the matrix.

## Usage

``` r
param_readable(s, eta, ...)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- ...:

  Passed to the method. No method in this package reads it.

## Value

`NULL` when the family declares nothing, otherwise a list with

- `value`:

  a named numeric vector of the quantities, on their own interpretable
  scale;

- `jacobian`:

  a numeric matrix with one row per quantity and one column per free
  value, row names matching `value` and no column names;

- `transform`:

  a character vector, one entry per quantity and named like `value`,
  naming the scale its interval is built on: one of `"identity"`,
  `"log"`, `"atanh"` or `"logit"`;

- `label`:

  a single string naming the block, for a consumer laying out a printed
  summary. The family supplies it because the family is what holds the
  reading.

## What a consumer does with it

Given the free vector's variance matrix \\V\\, the standard errors of
the declared quantities are the delta method, \\\sqrt{\mathrm{diag}(J V
J^\top)}\\. Each interval is then built on the scale `transform` names
and mapped back, so a variance stays positive and a correlation stays
inside \\(-1, 1)\\: on the raw scale an interval for a correlation
routinely runs past 1.

## The Jacobians are closed form

For a family whose coordinates are separate scalar links of separate
free values, \\J\\ is the diagonal matrix of the inverse links' first
derivatives, which
[`readable_diagonal()`](https://statmodels7.github.io/parameters7/reference/readable_diagonal.md)
assembles. For
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
it is the softmax Jacobian. For
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
the autoregressive coefficients come out of the Levinson-Durbin
recursion the family already propagates derivative arrays through, so
their derivatives in every free value are read off the first-order block
instead of being computed again.

## Which families declare something

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
and
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
declare a variance and a correlation.
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
declares its multiplier, or `NULL` when it was built with `link = NULL`
and has none.
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
declares the probability vector.
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
declares the variance, the partial autocorrelations that parametrize it,
**and** the autoregressive coefficients, which appear nowhere in the
covariance. Every other family, including
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
takes the base method and returns `NULL`.

## See also

[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
for the matrix itself,
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for the family with the most to declare, and
[`readable_diagonal()`](https://statmodels7.github.io/parameters7/reference/readable_diagonal.md)
for the shared assembly.

## Examples

``` r
# An AR(1) is about a variance and a correlation, not about a log and an
# inverse hyperbolic tangent.
r <- param_readable(ar1(5), c(log(2), atanh(0.6)))
r$value
#> scale   rho 
#>   2.0   0.6 
r$transform
#>   scale     rho 
#>   "log" "atanh" 

# What a consumer does with the Jacobian: the delta method. With a variance
# of 0.01 on the second free value, the correlation's standard error is
# |drho/dz| * 0.1.
V <- diag(c(0.04, 0.01))
sqrt(diag(r$jacobian %*% V %*% t(r$jacobian)))
#> scale   rho 
#> 0.400 0.064 

# An autoregression declares the coefficients as well as the partial
# autocorrelations, and at q = 2 they differ.
a <- param_readable(autoregressive(8, 2), c(0, 0.9, -0.4))
a$value
#>      scale      pacf1      pacf2       phi1       phi2 
#>  1.0000000  0.7162979 -0.3799490  0.9884545 -0.3799490 

# At q = 1 they coincide: the first partial autocorrelation IS the AR(1)
# coefficient.
b <- param_readable(autoregressive(5, 1), c(0, 0.9))
b$value[c("pacf1", "phi1")]
#>     pacf1      phi1 
#> 0.7162979 0.7162979 

# A family whose matrix is all there is to say declares nothing.
param_readable(log_cholesky(3), rep(0, 6))
#> NULL
```
