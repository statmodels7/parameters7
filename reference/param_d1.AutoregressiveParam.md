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
exactly and nothing is differenced. Measured against one central
difference of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
the first derivatives agree to \\5 \times 10^{-11}\\, which is the
difference's own accuracy.

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
and differ only in the order they pass.
[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
fills all four orders in one pass whatever is asked, so the first order
costs nearly what the fourth costs; see
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for the table.

## See also

[`ar_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_derivative.md),
which assembles them, and
[`param_dlogdet.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md)
for the log-determinant's own derivatives.
