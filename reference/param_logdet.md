# Log-Determinant of a Parameter's Matrix

Returns \\\log\|M\|\\ when the family is of full rank, and the log
pseudo-determinant – the sum of the logs of the non-zero eigenvalues –
when it is not.

## Usage

``` r
param_logdet(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Passed to methods.

## Value

A single number.

## Details

One generic covers both because a consumer asks the same question of
either: what normalising constant does this matrix contribute. The
object knows which answer is the right one, through its declared rank.

The sign in front of the result is the consumer's arithmetic. A gaussian
log-density written in the covariance carries
\\-\frac{1}{2}\log\|\Sigma\|\\ and one written in the precision carries
\\+\frac{1}{2}\log\|\Omega\|\\; the parameter answers what the
log-determinant is, and the likelihood decides where it goes.

## See also

[`param_dlogdet`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md),
[`param_solve`](https://statmodels7.github.io/parameters7/reference/param_solve.md)

## Examples

``` r
param_logdet(log_cholesky(2), c(0, 0, 0))
#> [1] 0

# a rank-deficient precision: the pseudo-determinant, over the 4 non-zero
# eigenvalues of a second-difference penalty on 6 coefficients
p <- crossprod(diff(diag(6), differences = 2))
s <- scaled_matrix(p)
c(rank = s@rank, logdet = param_logdet(s, 0))
#>    rank  logdet 
#> 4.00000 4.65396 
```
