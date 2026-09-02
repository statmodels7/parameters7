# Log-Determinant of an Inverse Parameter

\\\log\lvert S^{-1}\rvert = -\log\lvert S\rvert\\, so the value and
every derivative order are the inner family's negated. No determinant of
the inverted matrix is computed.

## Arguments

- s:

  An
  [`InverseParam()`](https://statmodels7.github.io/parameters7/reference/InverseParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
a single number; the derivative methods a named list or vector shaped as
the inner family's, with every entry negated.

## See also

[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
for the family.
