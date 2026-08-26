# Second Derivatives of a Correlation Parameter

Closed form, by the Leibniz rule on \\R = LL^\top\\ taken twice,

\$\$\partial\_{kl} R = L\_{kl} L^\top + L_k L_l^\top + L_l L_k^\top + L
L\_{kl}^\top,\$\$

with each \\\partial^S L\\ from
[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)
and zero wherever the multiset \\S\\ spans two rows of \\L\\.

What survives is worth knowing before reading a result. A component in
two angles from **different** rows \\i\\ and \\j\\ does not vanish: the
two outer terms drop, but \\L_k L_l^\top + L_l L_k^\top\\ is supported
on exactly the entries \\(i, j)\\ and \\(j, i)\\. It is zero even there
when the two angles sit beyond the columns the two rows share. The
diagonal is exactly zero at every order.

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 1, 2)` symmetric matrices keyed as
`param_tuple_names(s)` and in that order, each with a zero diagonal.

## See also

[`corr_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_derivative.md),
which assembles it,
[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)
for the vanishing rules, and
[`param_d1.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.CorrelationParam.md)
for the order below.
