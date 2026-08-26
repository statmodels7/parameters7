# Value of an Autoregressive Parameter

Returns the Toeplitz matrix \\\gamma_0 \rho\_{\lvert i-j \rvert}\\, the
autocorrelations coming from the Levinson-Durbin recursion of
[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
and the marginal variance from the scale link. Positive definite at
every free vector, the partial autocorrelations being inside \\(-1, 1)\\
by construction.

The cost is the recursion, which is linear in \\p\\: 0.00018 s at \\q =
1, p = 10\\ and 0.00143 s at \\q = 1, p = 200\\.

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

A positive definite Toeplitz `s@dimension` by `s@dimension` numeric
matrix with dimnames `v1`, `v2`, ...

## See also

[`param_free.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_free.AutoregressiveParam.md)
for the inverse, and
[`ar_assemble()`](https://statmodels7.github.io/parameters7/reference/ar_assemble.md),
which fills it.
