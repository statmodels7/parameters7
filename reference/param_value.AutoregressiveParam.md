# Value of an Autoregressive Parameter

Returns the Toeplitz matrix \\\gamma_0 \rho\_{\lvert i-j \rvert}\\, the
autocorrelations coming from the Levinson-Durbin recursion of
[`ar_tables()`](https://statmodels7.github.io/parameters7/reference/ar_tables.md)
and the marginal variance from the scale link. Positive definite at
every free vector in exact arithmetic, the partial autocorrelations
being inside \\(-1, 1)\\ by construction; see
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for the edge of the chart in double precision.

The cost is the recursion, which is linear in \\p\\, and the filling of
the matrix.

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
[`ar_tables()`](https://statmodels7.github.io/parameters7/reference/ar_tables.md),
which computes it.
