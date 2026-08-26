# Fourth Derivatives of a Log-Cholesky Parameter

Closed form, by the same Leibniz rule on \\M = L L^\top\\ as at third
order, with the four differentiations distributed over the two factors.
The derivative of \\L\\ in a multiset of indices is non-zero only for an
empty multiset, a single index, or a repetition of one diagonal free
value, so of the \\2^4\\ ways to split a quadruple between the factors
almost all contribute nothing.

Exact at this order, which is where the contract stops: a fourth-order
chain rule through a link needs this much. Nothing is differenced, so
the accuracy is machine precision, against the \\10^{-4}\\ a stencil
would give.

A component repeating a below-diagonal free value three or more times is
exactly zero, \\M\\ being quadratic in each of those.

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

A list of `choose(s@n_free + 3, 4)` symmetric matrices keyed as
`param_tuple_names(s, 4)` and in that order. At \\p = 3\\ that is 126
components, at \\p = 8\\ it is 82251.

## See also

[`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md),
which assembles it,
[`param_d3.LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.LogCholeskyParam.md)
for the order below, and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
for what a family without a closed form gets instead.
