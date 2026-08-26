# Log-Determinant Gradient of a Scaled Parameter

Closed form: \\r\\ h'(\eta)/h(\eta)\\, which under the default log link
is the **rank itself**, at every scale. That is the fact a penalty rests
on: it is the term that makes the smoothing parameter estimable, the
stationary point of \\\tfrac{\lambda}{2}\beta^\top P\beta -
\tfrac{r}{2}\log\lambda\\ being \\\lambda = r/(\beta^\top P \beta)\\.

The constant \\\log\|P\|\_+\\ contributes nothing, so a deficient family
answers with its rank, never with its dimension.

## Arguments

- s:

  A
  [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`, or empty
for a fixed parameter.

## See also

[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
for the estimability argument in full, and
[`param_d2logdet.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.ScaledMatrixParam.md)
for the order above.
