# Numerical Fourth Derivatives of a Parameter's Value

Estimates the distinct fourth derivatives, one product stencil per index
tuple, applied directly to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md).
It is the order-four analogue of
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
and the method
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
dispatches to when a family has not written its own. At this order a
difference is the least accurate of the four, which is why every family
in this package carries a closed form instead.

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
attainable accuracy of order \\\varepsilon^{1/3}\\. That is enough to
catch a transcription error in a closed form and not enough for use in a
fit.

## See also

[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md),
the generic that this serves,
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
for the order below, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
which uses this stencil as the reference for a family that is not a
matrix.

## Examples

``` r
# A scalar matrix: every order is the matrix again, so the error shows.
s <- scalar_matrix(2)
max(abs(numerical_d4(s, 0.3)[[1]] - param_value(s, 0.3)))
#> [1] 7.178321e-06

# Against a family that writes its own.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
ana <- param_d4(q, eta)
c(gap = max(abs(unlist(numerical_d4(q, eta)) - unlist(ana))),
  scale = max(abs(unlist(ana))))
#>          gap        scale 
#> 1.164773e-04 2.386920e+01 

# A 1 percent error in the largest entry would be about two thousand times
# the gap above, so a wrong closed form is caught.
0.01 * max(abs(unlist(ana)))
#> [1] 0.238692
```
