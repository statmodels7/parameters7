# The Spectral Decomposition of a Parameter's Matrix

Returns the eigenvalues and eigenvectors of \\M(\eta)\\, together with a
flag marking the directions that carry the matrix. The base-class
log-determinant and its first two derivatives read it, so one
decomposition serves them.

## Usage

``` r
param_spectrum(s, eta)
```

## Arguments

- s:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  object, whose `rank` and `dimension` are read.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list with four components

- `values`:

  numeric of length `s@dimension`, in decreasing order.

- `vectors`:

  a `s@dimension` by `s@dimension` orthonormal matrix, one eigenvector
  per column, matching `values`.

- `keep`:

  logical of length `s@dimension`, `TRUE` in the first `s@rank`
  positions.

- `matrix`:

  the value itself, as
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  returned it, so a caller needing both does not evaluate the map twice.

## Details

Which eigenvalues count as zero is settled **by position**, never by
size. [`eigen()`](https://rdrr.io/r/base/eigen.html) returns them in
decreasing order, so the first `s@rank` are kept and the rest dropped,
whatever their numerical values are. The rank itself comes from the
object and is never re-derived: counting eigenvalues above a relative
tolerance is not scale invariant, so a family whose components differ by
many orders of magnitude would be assigned a different rank at different
\\\eta\\, and a fitted model with smoothing parameters that far apart is
ordinary. The rank is fixed once, at construction, from the components;
see
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md).

The matrix is symmetrized as `(m + t(m)) / 2` before the decomposition,
so both triangles of an asymmetry of rounding size enter the
eigenvalues.

## See also

[`spectrum_pinv()`](https://statmodels7.github.io/parameters7/reference/spectrum_pinv.md),
which builds the Moore-Penrose inverse from this, and
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md),
which fixes the rank at construction.
