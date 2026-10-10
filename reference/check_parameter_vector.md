# The Reduced Battery for a Parameter That Is Not a Matrix

The checks that
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
runs for a family whose value is not a symmetric matrix and which
therefore has no log-determinant, no solve and no factor: the inverse
round trip and each of the four derivative orders against the
single-stencil numerical construction, and, for
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
that the value stays on the simplex and that every derivative component
sums to zero over the value index.

## Usage

``` r
check_parameter_vector(s, tol = 1e-06, verbose = TRUE)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object that is not a
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md).

- tol:

  The relative tolerance that a check has to meet. The thresholds of the
  four derivative orders are \\10^{-6}\\, \\10^{-5}\\, \\10^{-4}\\ and
  \\5 \times 10^{-3}\\ at the default `1e-6`, and each is multiplied by
  `tol / 1e-6` otherwise.

- verbose:

  Whether to print the table.

## Value

Invisibly, a data frame with columns `check`, `status` and `note`, with
seven rows for
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
and five for any other family. The `status` is `"OK"`, `"FAIL"` or, for
a derivative order that comes from a numerical fallback,
`"NOT CHECKED"`. The column `note` is **character**, holding a formatted
number, the string `numerical` or the empty string, where the matrix
branch of
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
returns a numeric `statistic`.

## Details

The last check is an identity the set itself supplies. Differentiating
\\\sum_a \pi_a = 1\\ gives \\\sum_a \partial \pi_a = 0\\, and
differentiating again gives the same for every higher order, so a
derivative array that fails the identity is wrong, whatever other check
it passes. A
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
satisfies it row by row, its rows being simplexes.

Three free vectors are used, drawn from `rnorm(s@n_free, sd = 0.8)`
after `set.seed(101)`, `set.seed(102)` and `set.seed(103)`, so the
results are exactly reproducible. The caller's own `.Random.seed` is
saved on entry and restored on exit through
[`capture_seed()`](https://statmodels7.github.io/parameters7/reference/capture_seed.md)
and
[`restore_seed()`](https://statmodels7.github.io/parameters7/reference/capture_seed.md),
so a call in the middle of a simulation leaves that simulation
unchanged. The derivative comparisons are against
[`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md)
through
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md),
whose accuracy decreases with the order, so the thresholds widen with
the order.

## See also

[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
which dispatches here, and
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
the two families in this package that reach it.
