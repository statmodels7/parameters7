# Free Vector of a Block Replication

Reads the free vector from the **first** diagonal block, by inverting
the inner parameter there, and then checks the whole matrix against the
replication that free vector implies. A matrix whose blocks are not
identical, or whose first block is outside the inner family, is
rejected.

The rejection is complete because the whole rebuild is compared, not the
blocks pairwise: the first block might invert cleanly while the fourth
is something else entirely.

## Arguments

- s:

  A
  [`KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/KronIdentityParam.md)
  object.

- m:

  A symmetric numeric matrix of side `s@dimension`, block diagonal with
  `s@param_params$m` identical blocks, already checked for shape and
  symmetry by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## See also

[`param_value.KronIdentityParam()`](https://statmodels7.github.io/parameters7/reference/param_value.KronIdentityParam.md),
the map this inverts.
