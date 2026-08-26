# Second Derivatives of a Diagonal Parameter

Closed form, from
[`linkfunctions7::d2linkinv()`](https://statmodels7.github.io/linkfunctions7/reference/d2linkinv.html).
The family is **separable**: each diagonal entry is a function of one
free value alone, so a component \\\partial\_{kl} M\\ with \\k \ne l\\
is exactly the zero matrix, and a pure one carries \\h''(\eta_k)\\ in
the entries that value owns. A
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
has one free value and so one component.

The same separability holds at third and fourth order, which is why
those methods keep only the pure components too.

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

A list of `choose(s@n_free + 1, 2)` diagonal matrices keyed as
`param_tuple_names(s)` and in that order, the mixed ones exactly zero.

## See also

[`param_d1.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d1.DiagMatrixParam.md)
and
[`param_d3.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.DiagMatrixParam.md)
for the neighboring orders.
