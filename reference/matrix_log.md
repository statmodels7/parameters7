# Construct a Matrix Logarithm Parameter

The matrix logarithm parametrisation of a symmetric positive definite
matrix: \\M = \exp(S)\\ with \\S\\ symmetric and genuinely free – no
entry is transformed, the free values fill the lower triangle of \\S\\
directly, diagonal first and then below the diagonal column by column.

## Usage

``` r
matrix_log(dimension, role = c("either", "covariance", "precision"))
```

## Arguments

- dimension:

  The side \\p\\ of the matrix.

- role:

  A label; see
  [`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Value

An object of class
[`MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md).

## Details

Two quantities are exact by construction and cost nothing. The
log-determinant is \\\log\|M\| = \mathrm{tr}(S)\\, the sum of the
diagonal free values, so it is linear and its higher derivatives vanish;
and the inverse is \\M^{-1} = \exp(-S)\\, through the same
eigendecomposition that evaluates the map, with no factorisation.

The derivatives of the map are the Frechet derivatives of the matrix
exponential, by the Daleckii-Krein representation: with \\S = Q \Lambda
Q^\top\\ and directions rotated by \\Q\\, the \\k\\-th derivative
contracts the directions against divided differences of \\e^x\\ of order
\\k+1\\ at the eigenvalues, summed over the orderings of the directions.
The divided differences are computed by the Opitz theorem – the
exponential of a small upper bidiagonal matrix, read off its corner –
which stays exact under repeated and nearly repeated eigenvalues, where
the quotient formula cancels catastrophically.

Next to
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md):
the log-Cholesky derivatives are sparse products and cheaper, while here
the log-determinant and the inverse are the free quantities. The matrix
logarithm is the chart to reach for when both the covariance and the
precision are wanted at once.

## References

Daleckii, J. L. and Krein, S. G. (1965). Integration and differentiation
of functions of Hermitian operators. *American Mathematical Society
Translations* 47, 1-30.

Opitz, G. (1964). Steigungsmatrizen. *Zeitschrift fur Angewandte
Mathematik und Mechanik* 44, T52-T54.

## See also

[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)

## Examples

``` r
s <- matrix_log(2)
eta <- c(0.2, -0.3, 0.4)
round(param_value(s, eta), 4)
#>        v1     v2
#> v1 1.3058 0.3948
#> v2 0.3948 0.8123

# the log-determinant is the trace of S, linear in the free vector
c(param_logdet(s, eta), 0.2 - 0.3)
#> [1] -0.1 -0.1

# the round trip closes exactly
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 1.665335e-16
```
