# The Reduced Battery for a Parameter That Is Not a Matrix

What
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
runs for a family whose value is not a symmetric matrix, so has no
log-determinant, no solve and no factor. Seven checks: the inverse round
trip, each of the four derivative orders against the single-stencil
numerical construction, that the value stays on the simplex, and that
every derivative component sums to zero over the value index.

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

  The relative tolerance a check has to meet.

- verbose:

  Whether to print the table.

## Value

Invisibly, a seven-row data frame with columns `check`, `status` and
`note`. Note that `note` is **character**, holding a formatted number or
the empty string, where
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)'s
matrix branch returns a numeric `statistic`.

## Details

The last check is an identity the set itself supplies. Differentiating
\\\sum_a \pi_a = 1\\ gives \\\sum_a \partial \pi_a = 0\\, and
differentiating again gives the same for every higher order, so a
derivative array that does not sum to zero is wrong whatever else it
agrees with. A
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
which at third and fourth order are themselves good to about
\\10^{-5}\\: measured on
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
the four orders come back at \\4 \times 10^{-12}\\, \\4 \times
10^{-12}\\, \\3 \times 10^{-7}\\ and \\3 \times 10^{-5}\\, and the
widening is the reference's.

## See also

[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
which dispatches here, and
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
the two families that reach it.
