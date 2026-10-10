# Solve of a Block Replication

Solves blockwise: \\M^{-1} = I_m \otimes S^{-1}\\, so the rows of `b`
belonging to each block are handed to the inner parameter's own
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
in turn. Nothing of side \\md\\ is factorized, and the method uses the
inner parameter's own route, including a closed inverse where it has
one.

The generic has already rejected a rank-deficient composite, which a
deficient inner parameter produces, and filled `b` with the identity
when the caller left it out.

## Arguments

- s:

  A
  [`KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/KronIdentityParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- b:

  A numeric matrix with `s@dimension` rows, already coerced from a
  vector and defaulted to the identity by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric matrix with `s@dimension` rows and as many columns as `b`.

## See also

[`param_factor.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_factor.KronIdentityParam.md)
for the lifted factor, and
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
for why the deficiency is inherited.
