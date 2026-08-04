# The Softmax Point Behind a Free Vector

The value \\\pi\\, computed through the log-sum-exp shift so that a
large free value saturates gracefully instead of overflowing.

## Usage

``` r
simplex_point(eta)
```

## Arguments

- eta:

  A numeric vector of length \\K - 1\\.

## Value

A numeric vector of length \\K\\.
