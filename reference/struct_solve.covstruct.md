# Default Solve

Fallback: a Cholesky of
[`struct_matrix`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md),
with the definiteness verdict taken spectrally (see
[`chol_pd`](https://statmodels7.github.io/covstructs7/reference/chol_pd.md)).

## Arguments

- s:

  A
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  object.

- eta:

  A numeric vector of free values.

- b:

  A numeric matrix with `s@dimension` rows.

- ...:

  Unused.

## Value

A numeric matrix.
