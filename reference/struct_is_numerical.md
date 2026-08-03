# Which of a Structure's Quantities Come From the Base Class

Reports, for each of the five derivative quantities, whether the
structure supplies its own method or falls back to the one registered on
[`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

## Usage

``` r
struct_is_numerical(s)
```

## Arguments

- s:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

## Value

A named logical vector over `struct_dmatrix`, `struct_d2matrix`,
`struct_logdet`, `struct_dlogdet` and `struct_d2logdet`, `TRUE` where
the base-class method is in force.

## Details

The distinction that matters is whether an independent check exists. A
derivative computed by finite differences cannot be checked against a
finite difference, and a log-determinant read off an eigendecomposition
cannot be checked against an eigendecomposition; the comparison is the
same arithmetic twice, and it agrees however wrong the structure is.
[`check_covstruct`](https://statmodels7.github.io/covstructs7/reference/check_covstruct.md)
uses this to report such a quantity as not checked rather than as
passed.

[`struct_solve()`](https://statmodels7.github.io/covstructs7/reference/struct_solve.md)
and
[`struct_factor()`](https://statmodels7.github.io/covstructs7/reference/struct_factor.md)
are deliberately absent. Their base-class versions are a Cholesky
factorisation, which is exact whoever performs it, and the validator
compares them with [`base::solve`](https://rdrr.io/r/base/solve.html)
either way. Calling them numerical would suggest an approximation that
is not there.

## Examples

``` r
# everything closed form except the second-order log-determinant of a
# structure that does not supply one
struct_is_numerical(log_cholesky(3))
#>  struct_dmatrix struct_d2matrix   struct_logdet  struct_dlogdet struct_d2logdet 
#>           FALSE           FALSE           FALSE           FALSE           FALSE 
struct_is_numerical(diag_struct(3))
#>  struct_dmatrix struct_d2matrix   struct_logdet  struct_dlogdet struct_d2logdet 
#>           FALSE           FALSE           FALSE           FALSE           FALSE 
```
