# Log-Determinant of a Sum of Fixed Matrices

The log-determinant of the assembled matrix where the family has full
rank, and its log **pseudo**-determinant, the sum over the `s@rank`
largest eigenvalues, where it does not. A sum of fixed matrices has no
structure a closed form could exploit, so this is the one family here
that assembles and decomposes.

## Arguments

- s:

  A
  [`SumStructParam()`](https://statmodels7.github.io/parameters7/reference/SumStructParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A single number.

## Details

Which branch is taken is decided by `s@rank`, fixed at construction from
the components' shared null space, and never by counting eigenvalues at
the point: that count falls as the weights spread apart, which
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
measures. So the number of eigenvalues summed does not move as a fit
walks the free vector.

Measured on `sum_struct(list(matrix(1, 3, 3)))`, of rank 1, at a weight
of 2: the pseudo-determinant is \\\log 6\\, its one non-zero eigenvalue,
to the printed digit.

## See also

[`param_dlogdet.SumStructParam()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.SumStructParam.md)
for its derivatives, and
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
for the directions the pseudo-determinant leaves out.
