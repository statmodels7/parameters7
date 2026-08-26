# The Shared Null Space of a Set of Positive Semidefinite Matrices

Returns an orthonormal basis of the **intersection** of the components'
null spaces, which is the null space of every non-negative combination
of them: a quadratic form \\v^\top M v = \sum_k c_k\\ v^\top P_k v\\ is
a sum of non-negative terms, so it vanishes only where every term does.
It is therefore a property of the family and can be computed once, at
construction.

## Usage

``` r
sum_struct_null_basis(components)
```

## Arguments

- components:

  A list of symmetric positive semidefinite matrices of the same side,
  already checked by
  [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md).

## Value

A numeric matrix with `nrow(components[[1]])` rows whose columns are an
orthonormal basis of the shared null space, with **zero columns** where
the family has full rank.

## Details

The components are stacked into one tall matrix and its right singular
vectors at negligible singular values are returned, which is the
intersection of the row spaces' orthogonal complements.

Each component is **divided by its own largest entry** before stacking.
Without that normalization a component whose scale is many orders below
another's sinks below the tolerance and is read as absent, which is
exactly the failure a rank taken from an assembled matrix shows: see
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
for the measured table, where the eigenvalue count of \\M(\eta)\\ falls
from 4 to 1 as the weights spread while the family's rank stays 4.

The tolerance is LAPACK's usual one, `max(dim) * eps * max(d)`. Where
the stacked matrix has fewer singular values than the side, the missing
directions are null by construction and are appended.

## See also

[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md),
the only caller, and
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
for what the basis is used for.
