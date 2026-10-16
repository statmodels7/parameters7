# Numerical Fourth Derivatives of a Parameter's Value

Estimates the distinct fourth derivatives, one tensor-product stencil
per index tuple, through
[`param_tensor()`](https://statmodels7.github.io/parameters7/reference/param_tensor.md),
applied to the highest analytic derivative below the fourth that the
family supplies and to
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
where there is none. It is the order-four analogue of
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
and the method
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
dispatches to when a family has not written its own.

## Usage

``` r
numerical_d4(s, eta, accuracy = 4L)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- accuracy:

  The order of accuracy of the central rules, a positive even whole
  number.

## Value

A list of `choose(s@n_free + 3, 4)` estimates keyed as
`param_tuple_names(s, 4)` and in that order, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## Details

The tensor product is the one written out under
[`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md).
With an analytic
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
each tuple is one first difference of a third derivative, and the result
is accurate to about \\10^{-12}\\. From
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
alone the difference is of order four: rounding is then of order
\\\varepsilon / h^4\\, and the attainable accuracy at the default
accuracy 4 is of order \\\varepsilon^{1/2}\\, enough to catch a
transcription error in a closed form.

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
#> [1] 1.063594e-13

# Against a family that writes its own.
q <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
ana <- param_d4(q, eta)
c(gap = max(abs(unlist(numerical_d4(q, eta)) - unlist(ana))),
  scale = max(abs(unlist(ana))))
#>          gap        scale 
#> 4.469314e-12 2.386920e+01 
```
