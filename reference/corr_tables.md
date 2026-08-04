# Sines and Cosines of a Correlation Parameter's Angles

For every angle, the value and the first four derivatives, in the free
value, of both \\\sin\theta\\ and \\\cos\theta\\.

## Usage

``` r
corr_tables(s, eta)
```

## Arguments

- s:

  A
  [`CorrelationParam`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A list with `sin` and `cos`, each a list indexed by free value holding
five numbers: the value and four derivatives.

## Details

Each angle depends on exactly one free value, so these tables are all
the derivative machinery the family needs: an entry of \\L\\ is a
product of such factors, and differentiating it replaces each factor by
the derivative of the matching order. The composition of the
trigonometric function with the link is done by
[`compose4`](https://statmodels7.github.io/parameters7/reference/compose4.md).
