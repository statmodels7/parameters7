# Numerical First Derivatives of a Structure's Matrix

One central difference of
[`struct_matrix`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md)
in each component of \\\eta\\.

## Usage

``` r
numerical_dmatrix(s, eta)
```

## Arguments

- s:

  A
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A named list of symmetric matrices.

## See also

[`struct_dmatrix`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.md)
