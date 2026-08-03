# Construct an Unstructured Positive Definite Structure

The log-Cholesky parametrisation of a symmetric positive definite
matrix: \\M = L L^\top\\ with \\L\\ lower triangular and positive on the
diagonal, the free values being the logarithms of the diagonal entries
of \\L\\ and the entries below it.

## Usage

``` r
log_cholesky(dimension, role = c("either", "covariance", "precision"))
```

## Arguments

- dimension:

  The side \\p\\ of the matrix.

- role:

  A label, one of `"either"` (the default), `"covariance"` or
  `"precision"`. Nothing computed depends on it; it records which side
  of a model the matrix parametrises, since the family name does not
  say.

## Value

An object of class
[`LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/LogCholeskyStruct.md).

## Details

The parametrisation is that of Pinheiro and Bates (1996), and it is the
one to reach for when nothing is known about the matrix. It is a
bijection onto the positive definite cone, smooth in both directions,
with no boundary to reach on the free scale, and unique because the
Cholesky factor with a positive diagonal is.

The free vector runs the logarithms of the diagonal first and then the
below-diagonal entries column by column, so for \\p = 3\\ it is \\(\log
L\_{11}, \log L\_{22}, \log L\_{33}, L\_{21}, L\_{31}, L\_{32})\\. The
ordering is fixed rather than incidental: `free_names` follows it and
every consumer builds its parameter tables from those names.

The log-determinant is linear in the free vector, \\\log\|M\| = 2\sum_i
\log L\_{ii}\\, so its gradient is the constant 2 in the diagonal
directions and 0 elsewhere, and its Hessian vanishes. Both are supplied
in closed form.

The log on the diagonal is intrinsic to the parametrisation and not a
swappable link, which is why it appears in the free names.

## References

Pinheiro, J. C. and Bates, D. M. (1996). Unconstrained parametrizations
for variance-covariance matrices. *Statistics and Computing* 6, 289-296.

## See also

[`diag_struct`](https://statmodels7.github.io/covstructs7/reference/diag_struct.md),
[`scaled_struct`](https://statmodels7.github.io/covstructs7/reference/scaled_struct.md),
[`check_covstruct`](https://statmodels7.github.io/covstructs7/reference/check_covstruct.md)

## Examples

``` r
s <- log_cholesky(3)
s
#> Structure: log_cholesky
#> Matrix:    3 x 3, symmetric
#> Rank:      3 of 3
#> Role:      either
#> 
#> Free values: 6
#>   log_L1, log_L2, log_L3, L2.1, L3.1, L3.2
#> 
#> From the base class: none

eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
round(struct_matrix(s, eta), 4)
#>         v1      v2      v3
#> v1  1.2214  0.5526 -0.4421
#> v2  0.5526  0.9203 -0.0363
#> v3 -0.4421 -0.0363  2.0221

# the round trip closes exactly
max(abs(struct_free(s, struct_matrix(s, eta)) - eta))
#> [1] 6.938894e-17

# and the log-determinant is linear in the free vector
struct_dlogdet(s, eta)
#> log_L1 log_L2 log_L3   L2.1   L3.1   L3.2 
#>      2      2      2      0      0      0 
```
