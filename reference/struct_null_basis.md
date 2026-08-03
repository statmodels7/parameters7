# An Orthonormal Basis of a Null Space, and the Rank That Goes With It

Returns the rank and an orthonormal basis of the common null space of a
set of symmetric positive semidefinite matrices.

## Usage

``` r
struct_null_basis(mats, tol = 1e-10)
```

## Arguments

- mats:

  A list of symmetric matrices, or a single matrix.

- tol:

  The relative tolerance below which a singular value counts as zero.

## Value

A list with `rank` and `null_basis`.

## Details

The components are normalised individually before being stacked, and the
rank is read off the stack rather than off any assembled combination.
The reason is that a rank read off a sum is not scale invariant: on a
tensor-product penalty of true rank 28 out of 32, counting eigenvalues
above a relative tolerance of \\\lambda_1 P_1 + \lambda_2 P_2\\ gives 28
when the two are equal, 26 at a ratio of \\10^8\\ and 16 at \\10^{12}\\,
because the smaller contributions sink below the tolerance and are read
as zeros. Smoothing parameters that far apart are an ordinary fitted
model.

The null space of a sum of positive semidefinite matrices is the
intersection of their null spaces, which is the null space of the stack,
so the stacked route is not an approximation of the assembled one but
the exact statement of the same quantity.

## See also

[`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)

## Examples

``` r
# a second-difference penalty: its null space is the constants and the lines
d <- diff(diag(6), differences = 2)
struct_null_basis(crossprod(d))$rank
#> [1] 4
```
