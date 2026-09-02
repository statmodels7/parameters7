# Construct an Unstructured Positive Definite Parameter

Estimating a covariance matrix directly is awkward. An optimizer moving
freely through \\p(p+1)/2\\ numbers will eventually propose a matrix
that is not positive definite, and the likelihood is undefined there.
This function removes the problem instead of policing it: it returns an
object holding the log-Cholesky map, which sends an unconstrained vector
to \\M = L L^\top\\ with \\L\\ lower triangular and positive on the
diagonal. Every vector in \\\mathbb{R}^{p(p+1)/2}\\ gives a valid
matrix, so the optimizer never has to be told about the constraint.

Reach for it when nothing is known about the matrix. It is the
parametrization of Pinheiro and Bates (1996).

## Usage

``` r
log_cholesky(dimension)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix. A single positive whole number, finite
  and at least 1. `0`, `2.5`, `c(1, 2)`, `"3"`, `Inf` and `NA` all throw
  `'dimension' must be a single positive integer.`

## Value

An object of class
[`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md),
with properties

- `dimension`:

  integer, the side \\p\\.

- `n_free`:

  integer, \\p(p+1)/2\\.

- `free_names`:

  character, `log_L1` ... `log_Lp` then `L2.1`, `L3.1`, `L3.2`, ... in
  the order above.

- `rank`:

  integer, always \\p\\: the map is onto the full-rank cone.

- `null_basis`:

  a \\p\\ by 0 matrix; there are no null directions.

- `param_name`:

  `"log_cholesky"`.

- `param_params`:

  a list with one entry, `positions`, holding the row, column and
  diagonal flag of each free value.

## Why the map is safe

\$\$M = L L^\top, \qquad L\_{ii} = \exp(\eta_i) \> 0\$\$

The exponential keeps the diagonal of \\L\\ positive, and a triangular
matrix with a positive diagonal has full rank, so \\L L^\top\\ is
positive definite at every finite \\\eta\\. The map is smooth in both
directions and one to one, the Cholesky factor with a positive diagonal
being unique, so there is no boundary on the free scale to run into.

The logarithm on the diagonal belongs to the parametrization and is not
a swappable link, which is why it appears in the free names. Use
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
if a choice of link is wanted.

## Reading the free vector

The \\p\\ logarithms of the diagonal come first, then the strictly
below-diagonal entries column by column. For \\p = 3\\:

\$\$\eta = (\log L\_{11}, \log L\_{22}, \log L\_{33}, L\_{21}, L\_{31},
L\_{32})\$\$

Its length is \\p(p+1)/2\\, which is 3, 6, 10 and 15 at \\p\\ of 2 to 5.
The ordering is part of the interface: `free_names` follows it and
downstream parameter tables are built from those names. Treat it as
fixed.

## The log-determinant comes free

A likelihood involving \\M\\ almost always needs \\\log\|M\|\\, and here
it is twice the sum of the first \\p\\ free values:

\$\$\log\|M\| = 2 \sum\_{i=1}^{p} \log L\_{ii} = 2 \sum\_{i=1}^{p}
\eta_i\$\$

Linear, so the gradient is the constant 2 in the diagonal directions and
0 elsewhere, and the second, third and fourth derivatives are exactly
zero. All four are returned without any factorization being taken.

## When to use something else

This is the parametrization for an unstructured matrix, where
\\p(p+1)/2\\ free values is the price of assuming nothing. A structured
family estimates fewer and imposes its structure exactly:
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
for independence,
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
for exchangeable equicorrelation,
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
and
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
for a time series,
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
for a correlation with a separate scale,
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
for the same cone through the matrix exponential.

## Notation

\\\eta\\ is the free vector, the point on the unconstrained scale that
an optimizer moves, of length \\d = p(p+1)/2\\. \\p\\ is the side of the
matrix, \\L\\ the lower triangular Cholesky factor and \\M = L L^\top\\
the matrix itself. \\E\_{ij}\\ is the matrix with 1 in position \\(i,
j)\\ and 0 elsewhere.

## References

Pinheiro, J. C. and Bates, D. M. (1996). Unconstrained parametrizations
for variance-covariance matrices. *Statistics and Computing* **6**,
289-296.

## See also

[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
for the map and its inverse,
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md),
which returns \\L\\ itself at no cost,
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
and
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
for the determinant, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
to verify a parametrization.

[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
and
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
for the alternatives.

## Examples

``` r
# A 3 x 3 unstructured covariance: six free values, no constraints.
s <- log_cholesky(3)
s
#> Parameter: log_cholesky
#> Matrix:    3 x 3, symmetric
#> Rank:      3 of 3
#> 
#> Free values: 6
#>   log_L1, log_L2, log_L3, L2.1, L3.1, L3.2
#> 
#> From the base class: none
s@n_free
#> [1] 6
s@free_names
#> [1] "log_L1" "log_L2" "log_L3" "L2.1"   "L3.1"   "L3.2"  

# Pick any six numbers at all; the result is positive definite.
eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
M <- param_value(s, eta)
round(M, 4)
#>         v1      v2      v3
#> v1  1.2214  0.5526 -0.4421
#> v2  0.5526  0.9203 -0.0363
#> v3 -0.4421 -0.0363  2.0221
eigen(M, only.values = TRUE)$values
#> [1] 2.2750769 1.4303176 0.4584471

# Even absurd values stay in the cone, which is the point.
eigen(param_value(s, c(-8, 9, -7, 100, -100, 50)),
      only.values = TRUE)$values > 0
#> [1] TRUE TRUE TRUE

# The inverse map recovers the free vector exactly.
max(abs(param_free(s, M) - eta))
#> [1] 6.938894e-17

# log|M| is linear in eta: twice the sum of the first p entries.
all.equal(param_logdet(s, eta), 2 * sum(eta[1:3]))
#> [1] TRUE
all.equal(param_logdet(s, eta),
          sum(log(eigen(M, only.values = TRUE)$values)))
#> [1] TRUE
param_dlogdet(s, eta)
#> log_L1 log_L2 log_L3   L2.1   L3.1   L3.2 
#>      2      2      2      0      0      0 

# So every higher derivative of the log-determinant vanishes.
c(second = max(abs(param_d2logdet(s, eta))),
  third = max(abs(param_d3logdet(s, eta))),
  fourth = max(abs(param_d4logdet(s, eta))))
#> second  third fourth 
#>      0      0      0 

# The factor is what the parametrization holds, so it costs nothing.
L <- param_factor(s, eta)
round(L, 4)
#>         [,1]   [,2]   [,3]
#> [1,]  1.1052 0.0000 0.0000
#> [2,]  0.5000 0.8187 0.0000
#> [3,] -0.4000 0.2000 1.3499
all.equal(diag(L), exp(eta[1:3]))
#> [1] TRUE
```
