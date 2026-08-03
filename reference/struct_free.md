# The Free Vector Behind a Matrix

The inverse map: given a matrix in the family's set, returns the free
vector that produces it.

## Usage

``` r
struct_free(s, m, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

- m:

  A symmetric numeric matrix of the structure's dimension.

- ...:

  Passed to methods.

## Value

A numeric vector of length `s@n_free`, named by `s@free_names`.

## Details

Exact or refused, never obtained by optimisation. A generic
optimisation-based inverse would return a plausible \\\eta\\ for a
matrix outside the set, which is a wrong answer wearing the shape of a
right one. A family whose map has no closed-form inverse refuses
instead, and the base class refuses on behalf of any structure that does
not implement this.

## See also

[`struct_matrix`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md)

## Examples

``` r
s <- log_cholesky(2)
eta <- c(0.3, -0.2, 0.5)
struct_free(s, struct_matrix(s, eta))
#> log_L1 log_L2   L2.1 
#>    0.3   -0.2    0.5 
```
