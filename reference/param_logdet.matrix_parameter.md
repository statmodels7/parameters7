# Default Log-Determinant

The method every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
of its own. It takes an eigendecomposition of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
keeps the first `s@rank` eigenvalues, and returns the sum of their
logarithms: the log-determinant for a full-rank family and the log
pseudo-determinant otherwise. Which eigenvalues are kept is decided by
position, from the declared rank, and never from their size.

Exact, not approximated, so
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reporting `TRUE` here means that the answer costs \\O(p^3)\\ and cannot
be checked against an eigendecomposition, not that it is inaccurate.
Measured against a closed form on a \\4 \times 4\\ AR(1) covariance, the
agreement is \\4 \times 10^{-15}\\.

## Arguments

- s:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  object, whose `rank` and `dimension` are read.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number.

## Details

A non-positive eigenvalue among the ones the rank keeps throws, naming
the family and the counts: the family has declared a rank it does not
have at this \\\eta\\, so [`log()`](https://rdrr.io/r/base/Log.html) of
a non-positive number would be the wrong thing to return. This is the
check that catches a
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
method whose matrix leaves the positive semidefinite cone.

## See also

[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
for the generic,
[`param_spectrum()`](https://statmodels7.github.io/parameters7/reference/param_spectrum.md)
for the decomposition, and
[`param_dlogdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.matrix_parameter.md)
for its derivative.
