# Third Derivatives of a Diagonal Parameter

Closed form. A diagonal entry depends on one free value through its
link, so the only surviving components are the pure ones, and each
carries \\h'''(\eta_k)\\ on the entries that value owns. Every mixed
component is the exact zero matrix.

\\h'''\\ comes from
[`linkfunctions7::d3linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d3linkinv.html),
so the accuracy is the link's; a family without a closed form here would
get a product stencil good to about six digits instead.

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

A list of `choose(s@n_free + 2, 3)` diagonal matrices keyed as
`param_tuple_names(s, 3)` and in that order.

## See also

[`diag_higher()`](https://statmodels7.github.io/parameters7/reference/diag_higher.md),
which assembles it, and
[`param_d4.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.DiagMatrixParam.md)
for the order above.
