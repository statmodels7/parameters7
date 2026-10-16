# Third Derivatives of a Log-Cholesky Parameter

Computes the third derivatives in closed form, by the Leibniz rule on
\\M = L L^\top\\. Each component distributes its three differentiations
over the two factors, and a factor differentiated more than once
survives only where the repetitions name the same diagonal free value,
because every derivative of \\e^{\eta_k}\\ is itself.

A component repeating a below-diagonal free value three times is exactly
zero, because \\M\\ is quadratic in it. A component repeating a diagonal
value three times is not zero, because that value enters through an
exponential.

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 2, 3)` symmetric matrices keyed as
`param_tuple_names(s, 3)` and in that order.

## See also

[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
which assembles it, and
[`param_d4.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.LogCholeskyParam.md)
for the order above.
