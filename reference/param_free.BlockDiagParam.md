# Free Vector of a Block-Diagonal Parameter

Inverts each diagonal block through its own parameter and concatenates
the results. A matrix whose off-diagonal blocks are not zero is
rejected, that being a matrix the family cannot represent.

## Arguments

- s:

  A
  [`BlockDiagParam`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md)
  object.

- m:

  A symmetric matrix of side `s@dimension`.

- ...:

  Unused.

## Value

A numeric vector of length `s@n_free`.
