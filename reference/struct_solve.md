# Solve Through a Structure's Matrix

Returns \\M^{-1} B\\, computed through a factor rather than through an
explicit inverse.

## Usage

``` r
struct_solve(s, eta, b = NULL, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

- eta:

  A numeric vector of length `s@n_free`.

- b:

  A numeric matrix or vector with `s@dimension` rows. Defaults to the
  identity, which returns the inverse.

- ...:

  Passed to methods.

## Value

A numeric matrix with `s@dimension` rows.

## Details

A rank-deficient structure refuses rather than returning a
pseudo-inverse. What a consumer of an improper prior needs is the
quadratic form and the log pseudo-determinant – the penalised normal
equations invert \\X^\top X + \lambda P\\, which is non-singular even
when \\P\\ is not, and is assembled by the consumer – so a
pseudo-inverse would be a plausible matrix answering a question nobody
asked.

## See also

[`struct_factor`](https://statmodels7.github.io/covstructs7/reference/struct_factor.md),
[`struct_logdet`](https://statmodels7.github.io/covstructs7/reference/struct_logdet.md)

## Examples

``` r
s <- log_cholesky(2)
round(struct_solve(s, c(0, 0, 0.5)), 4)
#>       [,1] [,2]
#> [1,]  1.25 -0.5
#> [2,] -0.50  1.0
```
