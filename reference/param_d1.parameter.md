# Default First Derivatives

The method that every
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
inherits when it registers no
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
of its own. It estimates \\\partial V/\partial \eta_k\\ by one central
difference of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
in each component, through
[`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md)
at its default accuracy 4: the five-point rule at the step
\\\varepsilon^{1/5}\max(1, \|\eta_k\|)\\, at a cost of \\4d\\
evaluations of the map. The families in this package do not reach it.

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
the result of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and symmetrized for a matrix family.

## See also

[`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md),
which does the work, and
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
to find out whether a family reaches this method.
