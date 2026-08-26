# First Derivatives of a Matrix Logarithm Parameter

The Frechet derivative of the matrix exponential, by Daleckii-Krein:
rotate the basis direction \\E_k\\ by the eigenvectors of \\S\\, weight
it entrywise by the two-point divided differences \\e\[\lambda_i,
\lambda_j\]\\, and rotate back. Exact, with no differencing.

Note that \\\partial_k M\\ is **not** \\E_k \exp(S)\\: the exponential
of a matrix does not commute with an arbitrary direction, which is why
the weighting by divided differences appears at all. Where the
eigenvalues are equal the weight is \\e^{\lambda}\\ and the formula
reduces to the scalar one, which
[`dd_exp()`](https://statmodels7.github.io/parameters7/reference/dd_exp.md)'s
Opitz route delivers exactly.

Measured against a central difference of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
at \\p = 3\\, the agreement is \\1.5 \times 10^{-10}\\, which is the
difference's own accuracy.

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `s@n_free` symmetric matrices named by `s@free_names`, each
`s@dimension` by `s@dimension`.

## See also

[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md),
which does the contraction,
[`dd_exp()`](https://statmodels7.github.io/parameters7/reference/dd_exp.md)
for the divided differences, and
[`param_d2.MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/param_d2.MatrixLogParam.md)
for the order above.
