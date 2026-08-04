# Derivative Components of an Economical Parameter

Assembles one derivative order of \\M = \sigma^2 P(\rho)\\ from the
scale's derivatives and the pattern's.

## Usage

``` r
econ_derivative(s, eta, order, pattern)
```

## Arguments

- s:

  The parameter.

- eta:

  A numeric vector of two free values.

- order:

  The derivative order, 1 to 4.

- pattern:

  A function of the parameter and the scalars, returning the pattern and
  its four derivatives in the second free value.

## Value

A named list of symmetric matrices.

## Details

The value is a product of a function of the first free value and a
function of the second, so a component with \\a\\ scale indices and
\\b\\ correlation indices is the \\a\\-th derivative of the scale times
the \\b\\-th derivative of the pattern. Nothing is approximated and no
order is special.
