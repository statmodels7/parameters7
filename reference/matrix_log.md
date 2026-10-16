# Construct a Matrix Logarithm Parameter

Returns an object holding the matrix logarithm parametrization of a
symmetric positive definite matrix: \\M = \exp(S)\\ with \\S\\ symmetric
and **genuinely free**. No entry is transformed; the free values fill
the lower triangle of \\S\\ directly, the diagonal first and then below
the diagonal column by column, and the exponential does the rest. Any
vector in \\\mathbb{R}^{p(p+1)/2}\\ gives a positive definite matrix,
the exponential of a symmetric matrix having positive eigenvalues, as
long as no eigenvalue of \\S\\ exceeds about 709.78, above which the
exponential overflows in double precision.

It parametrizes the same cone as
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
with the same number of free values. The choice between the two depends
on which quantities should have a closed form; see **Details**.

## Usage

``` r
matrix_log(dimension)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix. A single positive whole number, finite
  and at least 1; any other value signals the error
  `'dimension' must be a single positive integer.`

## Value

An object of class
[`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md),
with `n_free` equal to \\p(p+1)/2\\, `free_names` `S1` ... `Sp` then
`S2.1`, `S3.1`, ..., `rank` equal to `dimension`, an empty `null_basis`,
`param_name` `"matrix_log"`, and `param_params` holding `positions`.

## The log-determinant and the inverse

\$\$\log\|M\| = \mathrm{tr}(S) = \sum\_{i=1}^{p} \eta_i,\$\$

the sum of the diagonal free values, so the log-determinant is linear
and its second, third and fourth derivatives are exactly zero. And the
inverse is \\M^{-1} = \exp(-S)\\, evaluated through the same
eigendecomposition, with no factorization, so `param_solve(s, eta)` and
`param_value(s, -eta)` agree up to rounding.

## The derivatives, and the reason for the Opitz route

They are the Frechet derivatives of the matrix exponential, by the
Daleckii-Krein representation. With \\S = Q \Lambda Q^\top\\ and the
directions rotated by \\Q\\, the \\k\\-th derivative contracts the
directions against divided differences of \\e^x\\ of order \\k+1\\ at
the eigenvalues, summed over the orderings of the directions.

The divided differences come from the Opitz theorem, the exponential of
a small upper bidiagonal matrix read off its corner, in place of the
recursive quotient, which cancels catastrophically when two eigenvalues
nearly coincide. The Opitz value stays accurate to rounding there.
Repeated eigenvalues occur in ordinary use: an \\S\\ proportional to the
identity, or any \\S\\ with a symmetry, has them exactly.

## Choosing between the two unstructured charts

The derivatives of
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
are sparse products of single-entry matrices and are far cheaper. Here a
fourth-order component sums over up to \\4! = 24\\ orderings of its
directions, each a contraction of \\O(p^5)\\ operations.

This chart has a linear log-determinant and an inverse that is the map
at \\-\eta\\, while the inverse of
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
needs a triangular solve. It suits a model that needs the covariance and
the precision together, or one whose free values should form an
unconstrained symmetric matrix.
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
suits a model that takes derivatives repeatedly.

## Notation

\\\eta\\ is the free vector, of length \\d = p(p+1)/2\\, \\S\\ the
symmetric matrix it fills, and \\M = \exp(S)\\ the value. \\Q\\ and
\\\Lambda\\ are the eigenvectors and eigenvalues of \\S\\, and
\\e\[\lambda_1, \dots, \lambda_m\]\\ a divided difference of the
exponential.

## References

Daleckii, J. L. and Krein, S. G. (1965). Integration and differentiation
of functions of Hermitian operators. *American Mathematical Society
Translations* **47**, 1-30.

Opitz, G. (1964). Steigungsmatrizen. *Zeitschrift fur Angewandte
Mathematik und Mechanik* **44**, T52-T54.

## See also

[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
for the other unstructured chart onto the same cone,
[`dd_exp()`](https://statmodels7.github.io/parameters7/reference/dd_exp.md)
for the divided differences, and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md),
which here is the map at \\-\eta\\.

## Examples

``` r
# A 3 x 3 unstructured covariance. The free names carry no transformation:
# the free values are the entries of S itself.
s <- matrix_log(3)
s@free_names
#> [1] "S1"   "S2"   "S3"   "S2.1" "S3.1" "S3.2"
eta <- c(0.2, -0.3, 0.4, 0.1, -0.2, 0.15)
M <- param_value(s, eta)
round(M, 4)
#>         v1     v2      v3
#> v1  1.2518 0.0805 -0.2650
#> v2  0.0805 0.7550  0.1517
#> v3 -0.2650 0.1517  1.5323
eigen(M, only.values = TRUE)$values > 0
#> [1] TRUE TRUE TRUE

# The log-determinant is the trace of S, so it is the sum of the first p
# free values, and it agrees with the eigenvalues.
c(closed = param_logdet(s, eta), trace = sum(eta[1:3]),
  from_eigen = sum(log(eigen(M, only.values = TRUE)$values)))
#>     closed      trace from_eigen 
#>        0.3        0.3        0.3 

# Which makes its gradient a vector of ones and zeros, and every higher
# order exactly zero.
param_dlogdet(s, eta)
#>   S1   S2   S3 S2.1 S3.1 S3.2 
#>    1    1    1    0    0    0 
c(second = max(abs(param_d2logdet(s, eta))),
  third = max(abs(param_d3logdet(s, eta))),
  fourth = max(abs(param_d4logdet(s, eta))))
#> second  third fourth 
#>      0      0      0 

# The inverse is the map at -eta, computed with no factorization.
max(abs(param_solve(s, eta) - param_value(s, -eta)))
#> [1] 1.554312e-15
max(abs(param_solve(s, eta) - solve(M)))
#> [1] 1.110223e-15

# The round trip closes.
max(abs(param_free(s, M) - eta))
#> [1] 1.831868e-15
```
