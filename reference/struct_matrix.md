# The Matrix a Structure Produces

Maps the free vector \\\eta\\ to the matrix it parametrises.

## Usage

``` r
struct_matrix(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md).

- eta:

  A numeric vector of length `s@n_free`.

- ...:

  Passed to methods.

## Value

A symmetric numeric matrix with `s@dimension` rows and columns.

## Details

This is the only generic a structure must implement. Everything else in
the package has a numerical method registered on the
[`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
class and is therefore available from this one alone.

The generic validates the free vector before dispatching, so every
method, including one written outside the package, refuses a vector of
the wrong length or one containing a non-finite value. The free scale
has no boundary to reach, so a non-finite entry is a defect in the
caller rather than a point of the domain.

## See also

[`struct_free`](https://statmodels7.github.io/covstructs7/reference/struct_free.md),
[`struct_dmatrix`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.md),
[`struct_logdet`](https://statmodels7.github.io/covstructs7/reference/struct_logdet.md)

## Examples

``` r
s <- log_cholesky(2)
round(struct_matrix(s, c(0, 0, 0.5)), 4)
#>     v1   v2
#> v1 1.0 0.50
#> v2 0.5 1.25
```
