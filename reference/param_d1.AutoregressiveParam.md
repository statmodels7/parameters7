# Derivatives of an Autoregressive Parameter

[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for an
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
parameter, closed form at every order. The map from the partial
autocorrelations to the matrix is polynomial, so the derivative arrays
propagated through the Levinson-Durbin recursion give each derivative
exactly and nothing is differenced.

Every component is Toeplitz, the structure being a property of the
family and fixed as the point moves.

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of symmetric `s@dimension` by `s@dimension` matrices with
dimnames, keyed as `param_tuple_names(s, order)` and in that order:
`s@n_free` of them at order 1, then `choose(s@n_free + k - 1, k)` at
order \\k\\.

## Details

The four share
[`ar_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_derivative.md)
and differ only in the order they pass. Each order runs its own kernel,
which returns the components of that order only.

## See also

[`ar_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_derivative.md),
which assembles them, and
[`param_dlogdet.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md)
for the log-determinant's own derivatives.
