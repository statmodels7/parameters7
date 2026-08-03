# Print a Covariance Structure

Prints the family, the shape of the matrix, the rank, the free values
and which of the derived quantities are numerical.

## Arguments

- x:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

- ...:

  Unused.

## Value

Invisibly `x`.

## Examples

``` r
log_cholesky(3)
#> Structure: log_cholesky
#> Matrix:    3 x 3, symmetric
#> Rank:      3 of 3
#> Role:      either
#> 
#> Free values: 6
#>   log_L1, log_L2, log_L3, L2.1, L3.1, L3.2
#> 
#> From the base class: none
scaled_struct(crossprod(diff(diag(6), differences = 2)))
#> Structure: scaled
#> Matrix:    6 x 6, symmetric
#> Rank:      4 of 6 (null space of dimension 2)
#> Role:      precision
#> 
#> Free values: 1
#>   scale
#> 
#> From the base class: none
```
