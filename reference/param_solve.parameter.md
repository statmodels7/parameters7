# Default Solve

Fallback: a Cholesky of
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md),
with the definiteness verdict taken spectrally (see
[`chol_pd`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)).

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  A numeric vector of free values.

- b:

  A numeric matrix with `s@dimension` rows.

- ...:

  Unused.

## Value

A numeric matrix.
