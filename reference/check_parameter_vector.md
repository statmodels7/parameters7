# The Reduced Battery for a Parameter That Is Not a Matrix

What
[`check_parameter`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
runs for a family whose value is not a symmetric matrix: the inverse
round trip, every derivative order against the single-stencil numerical
construction, and the identities the value's own set supplies – on the
simplex the entries sum to one and every derivative component sums to
zero over the value index, since differentiating \\\sum_a \pi_a = 1\\
kills every order; a transition matrix satisfies the same row by row.

## Usage

``` r
check_parameter_vector(s, tol = 1e-06, verbose = TRUE)
```

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- tol:

  The relative tolerance.

- verbose:

  Whether to print the table.

## Value

A data frame with one row per check, invisibly when printed.
