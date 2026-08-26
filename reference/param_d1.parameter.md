# Default First Derivatives

The method every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
of its own. It estimates \\\partial V/\partial \eta_k\\ by one
three-point central difference of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
in each component, \\(V(\eta + h e_k) - V(\eta - h e_k)) / 2h\\, at the
step \\\varepsilon^{1/3}\max(1, \|\eta_k\|)\\, which is about \\6.1
\times 10^{-6}\\ near the origin. It costs \\2d\\ evaluations of the map
and delivers about eleven digits; measured against a closed form on a
\\2 \times 2\\ covariance the gap is \\7 \times 10^{-11}\\ on entries of
size 3. No family in this package reaches it.

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, of any branch.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `s@n_free` estimates named by `s@free_names`, each shaped like
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)'s
result and symmetrized for a matrix family.

## See also

[`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md),
which does the work and carries the stencil, and
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
to ask whether a family reaches this method.
