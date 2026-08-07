# Quantities a Family Is About

The interpretable quantities of a parameter, with the Jacobian of the
map from the free vector and the scale each quantity's interval belongs
on.

## Usage

``` r
param_readable(s, eta, ...)
```

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values.

- ...:

  Passed to methods.

## Value

`NULL` when the family declares nothing, otherwise a list with

- `value`:

  a named numeric vector of the quantities;

- `jacobian`:

  a matrix with one row per quantity and one column per free value;

- `transform`:

  a character vector naming the scale each interval is built on, one of
  `"identity"`, `"log"`, `"atanh"` or `"logit"`;

- `label`:

  a single string naming the group, for a consumer laying out a printed
  summary. The family supplies it because the family is what holds the
  reading.

## Details

A consumer that has estimated the free vector and its variance matrix
\\V\\ obtains standard errors for these quantities by the delta method,
\\J V J'\\, and builds each interval on the declared scale before
mapping it back, so that a variance stays positive and a correlation
stays inside its interval. The Jacobians here are closed form: for the
families whose coordinates are separate scalar links it is the diagonal
of the inverse link's derivative, and for
[`autoregressive`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
the autoregressive coefficients carry their own derivatives out of the
Levinson-Durbin recursion, which already carries their derivatives.

The base class declares nothing, so a family that has no quantity beyond
the matrix it builds loses nothing by saying so: the matrix itself is
what a consumer already reports.

## See also

[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md),
[`autoregressive`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)

## Examples

``` r
param_readable(ar1(5), c(log(2), atanh(0.6)))$value
#> scale   rho 
#>   2.0   0.6 
param_readable(autoregressive(8, 2), c(0, 0.9, -0.4))$value
#>      scale      pacf1      pacf2       phi1       phi2 
#>  1.0000000  0.7162979 -0.3799490  0.9884545 -0.3799490 
```
