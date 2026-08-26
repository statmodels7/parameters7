# First Derivatives of a Parameter's Value

Differentiates the map once in each free value, returning \\\partial V /
\partial \eta_k\\ for \\k = 1, \dots, d\\ as a list of \\d\\ objects
each shaped like the value itself. A likelihood needs these to turn a
score in the matrix into a score in the coordinates an optimizer moves,
by the chain rule \\\partial \ell / \partial \eta_k = \sum\_{ij}
(\partial \ell / \partial M\_{ij}) (\partial M\_{ij} / \partial
\eta_k)\\.

## Usage

``` r
param_d1(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- eta:

  A numeric vector of length `s@n_free`, finite in every entry. Checked
  before dispatch; see
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  for the three rejections.

- ...:

  Passed to the method. No method in this package reads it.

## Value

A list of `s@n_free` derivatives, named by `s@free_names`. Each entry
has the shape of the value: a symmetric `s@dimension` by `s@dimension`
matrix for a matrix family, a numeric vector of length `s@n_free + 1`
for
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md).
The derivative of a symmetric matrix is symmetric, so each entry is too.

## The list is keyed by the free names

Entry `k` is the derivative in the free value `s@free_names[k]`, and the
list carries those names, so a consumer assembling a gradient can match
a coordinate to its label without tracking positions.

## Exact, or one stencil

Every family here writes this out in closed form. A family that does not
gets the method registered on
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md),
which applies **one** three-point central difference to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
in each component, at the step
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
gives for a first derivative. It is never a difference of a lower-order
difference, so the errors of two stages cannot compound.
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reports which route a given object takes, and returns `FALSE` throughout
for every family this package ships.

## Notation

\\\eta\\ is the free vector, of length \\d = \\ `s@n_free`, and \\V\\
the value
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns.

## See also

[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for the higher orders,
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
for the derivative of the log-determinant, and
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
to ask which route an object takes.

## Examples

``` r
s <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
d1 <- param_d1(s, eta)
names(d1)
#> [1] "log_L1" "log_L2" "L2.1"  
d1[["L2.1"]]
#>          v1       v2
#> v1 0.000000 1.221403
#> v2 1.221403 0.800000

# Each entry is a derivative, and a central difference of the map agrees.
h <- 1e-5
fd <- (param_value(s, eta + c(0, 0, h)) - param_value(s, eta - c(0, 0, h))) /
  (2 * h)
max(abs(d1[["L2.1"]] - fd))
#> [1] 1.235234e-12

# For a scaled precision the derivative is the matrix itself, the free
# value being the log of the scale.
P <- crossprod(diff(diag(4)))
r <- scaled_matrix(P)
all.equal(param_d1(r, 1.3)[[1]], param_value(r, 1.3))
#> [1] TRUE
```
