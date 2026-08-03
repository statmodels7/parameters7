# Construct a Scalar Multiple of the Identity

A diagonal matrix with one free value shared by every entry, so \\M =
\tau I\\ with \\\tau\\ positive.

## Usage

``` r
scalar_struct(
  dimension,
  link = linkfunctions7::log_link(),
  role = c("either", "covariance", "precision")
)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix.

- link:

  A linkfunctions7 link mapping the free scale to the positive scale.
  Defaults to
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html).

- role:

  A label; see
  [`log_cholesky`](https://statmodels7.github.io/covstructs7/reference/log_cholesky.md).

## Value

An object of class
[`DiagStruct`](https://statmodels7.github.io/covstructs7/reference/DiagStruct.md).

## Details

The simplest structure there is, and the one a random effect with a
single variance component uses. With the default log link it is the same
matrix as `scaled_struct(diag(dimension))`, reached from the other
direction.

## See also

[`diag_struct`](https://statmodels7.github.io/covstructs7/reference/diag_struct.md),
[`scaled_struct`](https://statmodels7.github.io/covstructs7/reference/scaled_struct.md)

## Examples

``` r
s <- scalar_struct(3)
c(n_free = s@n_free)
#> n_free 
#>      1 
round(struct_matrix(s, log(2)), 4)
#>    v1 v2 v3
#> v1  2  0  0
#> v2  0  2  0
#> v3  0  0  2
```
