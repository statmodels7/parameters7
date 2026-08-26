# Free Vector of an Autoregressive Parameter

Reads the marginal variance off the common diagonal entry and the
partial autocorrelations off the Levinson-Durbin recursion run forwards
on the autocorrelations. Exact where `m` is in the set: measured at \\p
= 6\\, \\q = 2\\, the round trip closes to \\1.1 \times 10^{-16}\\.

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- m:

  The covariance matrix of a stationary autoregression of order
  `s@param_params$order`, already checked for shape and symmetry by the
  generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

Four things are rejected rather than approximated, each with a message
of its own, because the set this family describes is much smaller than
the set of symmetric positive definite matrices and a caller who lands
outside it has a model error, not a rounding one:

- a diagonal that is not constant, which is not stationary;

- a matrix that is not Toeplitz, which is rejected instead of being
  averaged along its diagonals;

- an implied partial autocorrelation at or beyond \\\pm 1\\, or a
  vanishing prediction variance, neither of which is stationary;

- a Toeplitz matrix whose autocorrelations stop following the
  Yule-Walker recursion beyond lag \\q\\. This last one is what
  separates an AR(2) from an AR(3) at \\p \ge 4\\, and it is checked by
  rebuilding the matrix from the recovered free vector and comparing.

The tolerances are relative to `max(1, max(abs(m)))`: \\10^{-8}\\ for
the structure and \\10^{-7}\\ for the rebuild.

## See also

[`param_value.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_value.AutoregressiveParam.md),
the map this inverts.
