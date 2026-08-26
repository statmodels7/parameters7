# Third Derivatives of a Matrix Logarithm Parameter

Chains of three rotated directions contracted against four-point divided
differences, summed over the six orderings. Exact, where a numerical
route at this order keeps about six digits.

The cost grows quickly: six orderings per component, each an \\O(p^4)\\
contraction, and \\\binom{d+2}{3}\\ components. See
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
for the comparison with
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
whose derivatives are sparse products.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 2, 3)` symmetric matrices keyed as
`param_tuple_names(s, 3)` and in that order.

## See also

[`mlog_higher()`](https://statmodels7.github.io/parameters7/reference/mlog_higher.md),
which assembles it, and
[`param_d4.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.MatrixLogParam.md)
for the order above.
