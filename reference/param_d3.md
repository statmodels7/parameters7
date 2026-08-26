# Third Derivatives of a Parameter's Value

Differentiates the map three times, returning the distinct components
\\\partial^3 V / \partial \eta_k \partial \eta_l \partial \eta_m\\, each
shaped like the value. The derivative is symmetric in its indices, so
the \\\binom{d+2}{3}\\ unordered triples are returned, not the \\d^3\\
ordered ones. A third-order chain rule through a link needs them, and so
does the exact gradient of a marginal criterion, which differentiates a
penalized mode with respect to a hyperparameter.

## Usage

``` r
param_d3(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- eta:

  A numeric vector of length `s@n_free`, finite in every entry.

- ...:

  Passed to the method. No method in this package reads it.

## Value

A named list of `choose(s@n_free + 2, 3)` entries, keyed as
`param_tuple_names(s, 3)` and in that order. Each entry has the shape of
the value.

## Keys

The list is keyed by `param_tuple_names(s, 3)`, the lexicographic
combinations with repetition, and `param_tuple_indices(s, 3)` gives the
index triples in the same order. At \\d = 3\\ there are ten of them, at
\\d = 6\\ fifty-six.

## Exact, or one stencil

Every family here is closed form. The fallback on
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
applies one product stencil per triple, of the width each repeated index
calls for: an index appearing three times takes the four-point
third-difference factor, one appearing twice the three-point
second-difference factor, and a distinct index a two-point
first-difference factor. The product is evaluated in one pass over
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
so no order is reached by differencing the order below.

## Notation

\\\eta\\ is the free vector, of length \\d = \\ `s@n_free`, and \\V\\
the value
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns.

## See also

[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for the neighboring orders,
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
for the keys, and
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
for the third derivative of the log-determinant.

## Examples

``` r
# A scalar matrix is exp(eta) times the identity, so every order in the one
# free value is the matrix again.
s <- scalar_matrix(2)
param_d3(s, 0.3)
#> $`log_scale:log_scale:log_scale`
#>          v1       v2
#> v1 1.349859 0.000000
#> v2 0.000000 1.349859
#> 
all.equal(param_d3(s, 0.3)[[1]], param_value(s, 0.3))
#> [1] TRUE

# log_cholesky is quadratic in each below-diagonal free value, so a triple
# that repeats one of those three times vanishes exactly.
q <- log_cholesky(2)
d3 <- param_d3(q, c(0.2, -0.1, 0.4))
length(d3)
#> [1] 10
max(abs(d3[["L2.1:L2.1:L2.1"]]))
#> [1] 0

# It does not vanish in a diagonal free value, which enters through a log.
d3[["log_L1:log_L1:log_L1"]]
#>            v1        v2
#> v1 11.9345976 0.4885611
#> v2  0.4885611 0.0000000
```
