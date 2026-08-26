# Free Vector of a Sum of Fixed Matrices

Recovers the weights by least squares on the components' entries,
stacking the \\P_k\\ as columns and solving against the entries of `m`.
Exact where `m` is in the span: measured on two variance components over
three coefficients, the round trip closes to \\3 \times 10^{-16}\\.

## Arguments

- s:

  A
  [`SumStructParam()`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- m:

  A symmetric `s@dimension` by `s@dimension` matrix lying in the
  non-negative span of the components, already checked for shape and
  symmetry by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, the linked weights.

## Details

Two refusals, and they are different failures:

- a matrix the combination cannot reproduce, checked by residual at
  \\10^{-8}\\ relative to `max(1, max(abs(m)))`. The span is
  \\K\\-dimensional inside the \\p(p+1)/2\\-dimensional space of
  symmetric matrices, so almost every matrix is outside it and this is
  the ordinary failure.

- a matrix in the span needing a **non-positive weight**, which the link
  cannot carry. Such a matrix may still be positive semidefinite, so the
  refusal is about the parametrization, never about the matrix.

The least-squares solve is exact wherever it succeeds, the residual
check being what separates a solution from a projection.

## See also

[`param_value.SumStructParam()`](https://statmodels7.github.io/parameters7/reference/param_value.SumStructParam.md),
the map this inverts.
