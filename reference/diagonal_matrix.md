# Construct a Diagonal Parameter

A diagonal matrix whose entries are positive, each carried through a
scalar link from linkfunctions7.

## Usage

``` r
diagonal_matrix(
  dimension,
  link = linkfunctions7::log_link(),
  role = c("either", "covariance", "precision")
)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix.

- link:

  A linkfunctions7 link mapping the free scale to the positive entries.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- role:

  A label; see
  [`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Value

An object of class
[`DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md).

## Details

This is the case in which the two packages meet, and it is the reason
parameters7 composes with linkfunctions7 rather than competing with it.
The Jacobian of a diagonal block is diagonal, which is exactly the
contract a scalar link satisfies, so the link objects are reused as they
are and their exact derivatives come with them.

A single free value shared by every entry is
[`scalar_matrix`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md).

## See also

[`scalar_matrix`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
[`log_cholesky`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)

## Examples

``` r
s <- diagonal_matrix(3)
round(param_value(s, c(0, 0.5, -0.5)), 4)
#>    v1     v2     v3
#> v1  1 0.0000 0.0000
#> v2  0 1.6487 0.0000
#> v3  0 0.0000 0.6065

# any link whose range is positive works
round(param_value(diagonal_matrix(2, link = linkfunctions7::sqrt_link()),
                    c(1, 2)), 4)
#>    v1 v2
#> v1  1  0
#> v2  0  4
```
