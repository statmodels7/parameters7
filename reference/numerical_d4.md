# Numerical Fourth Derivatives of a Parameter's Value

Estimates the distinct fourth derivatives, one product stencil per index
tuple, applied directly to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md).
It is the order-four analogue of
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
and the method
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
dispatches to when a family has not written its own. Treat it as a
starting point: at this order a difference keeps about five digits,
which is why every family this package ships carries a closed form
instead.

## Usage

``` r
numerical_d4(s, eta)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list of `choose(s@n_free + 3, 4)` estimates keyed as
`param_tuple_names(s, 4)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## The stencil

The tensor product is the one written out under
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md),
with the multiplicities summing to four. A tuple naming one component
four times costs five evaluations of the map, one naming two components
twice each nine, and one naming four distinct components sixteen.

## Accuracy

Truncation is of order \\h^2\\ and rounding of order \\\varepsilon /
h^4\\, so balancing the two gives a step of \\\varepsilon^{1/6}\max(1,
\|\eta_k\|)\\, about \\2.5 \times 10^{-3}\\ near the origin, and an
attainable accuracy of order \\\varepsilon^{1/3}\\. Measured against the
closed form of a \\2 \times 2\\ log-Cholesky covariance, the agreement
is \\1.2 \times 10^{-4}\\ absolute on entries of size 24, so five
digits. That is enough to catch a transcription error in a closed form.
It is not enough to fit with.

## See also

[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md),
the generic this serves,
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
for the order below, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
which uses a stencil of this kind as the independent reference for a
closed form.

## Examples

``` r
# A scalar matrix: every order is the matrix again, so the error shows.
s <- scalar_matrix(2)
max(abs(numerical_d4(s, 0.3)[[1]] - param_value(s, 0.3)))
#> [1] 7.178321e-06

# Against a family that writes its own: five digits, as stated.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
ana <- param_d4(q, eta)
c(gap = max(abs(unlist(numerical_d4(q, eta)) - unlist(ana))),
  scale = max(abs(unlist(ana))))
#>          gap        scale 
#> 1.164773e-04 2.386920e+01 

# Which is still ample to catch a wrong closed form: a 1 per cent error in
# one component is four hundred times the noise.
0.01 * max(abs(unlist(ana)))
#> [1] 0.238692
```
