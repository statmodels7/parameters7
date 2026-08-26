# Fourth Derivatives of a Diagonal Parameter

Closed form, the same separable structure as at third order with
\\h''''(\eta_k)\\ in place of \\h'''(\eta_k)\\. Only the pure components
survive; every tuple naming two different free values is the exact zero
matrix.

\\h''''\\ comes from
[`linkfunctions7::d4linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d4linkinv.html).
This is the order at which a numerical route is least usable, keeping
about five digits, so the closed form matters most here.

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 3, 4)` diagonal matrices keyed as
`param_tuple_names(s, 4)` and in that order.

## See also

[`diag_higher()`](https://statmodels7.github.io/parameters7/reference/diag_higher.md),
which assembles it,
[`param_d3.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.DiagMatrixParam.md)
for the order below, and
[`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
for the alternative.
