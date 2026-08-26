# Print a Constrained Parameter

Prints a one-screen summary of a parametrization: the family name, the
shape and rank of the matrix, the role label, how many free values there
are and what they are called, and which of the nine derived quantities
come from the base class instead of a closed form.

## Arguments

- x:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

Invisibly `x`, so a print inside a pipe does not break it.

## Details

Three things about the output are worth knowing. The matrix block is
printed only for a
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
a
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
having no dimension or rank to report. The free names are truncated at
twelve, with a count of the rest, so a large unstructured covariance
does not fill the screen. And the last line, `From the base class`, is
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)'s
answer: `none` means every quantity is exact, and a list of generics
names the ones a stencil computes. Every family this package ships
prints `none`.

A rank-deficient family prints its null space's dimension beside the
rank, which is the quickest way to see that a penalty is improper.

## See also

[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md),
whose answer the last line reports, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
for a validation report instead of a description.

## Examples

``` r
# A full-rank family: nothing from the base class.
log_cholesky(3)
#> Parameter: log_cholesky
#> Matrix:    3 x 3, symmetric
#> Rank:      3 of 3
#> Role:      either
#> 
#> Free values: 6
#>   log_L1, log_L2, log_L3, L2.1, L3.1, L3.2
#> 
#> From the base class: none

# A rank-deficient one prints its null space beside the rank.
scaled_matrix(crossprod(diff(diag(6), differences = 2)))
#> Parameter: scaled
#> Matrix:    6 x 6, symmetric
#> Rank:      4 of 6 (null space of dimension 2)
#> Role:      precision
#> 
#> Free values: 1
#>   log_scale
#> 
#> From the base class: none

# A family that is not a matrix prints no matrix block.
simplex(4)
#> Parameter: simplex
#> 
#> Free values: 3
#>   alr1, alr2, alr3
#> 
#> From the base class: none

# Twelve free names at most, with the rest counted.
log_cholesky(6)
#> Parameter: log_cholesky
#> Matrix:    6 x 6, symmetric
#> Rank:      6 of 6
#> Role:      either
#> 
#> Free values: 21
#>   log_L1, log_L2, log_L3, log_L4, log_L5, log_L6, L2.1, L3.1, L4.1, L5.1, L6.1, L3.2, ... (9 more)
#> 
#> From the base class: none
```
