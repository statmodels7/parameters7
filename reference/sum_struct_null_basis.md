# The Shared Null Space of a Set of Positive Semidefinite Matrices

Returns an orthonormal basis of the intersection of the components' null
spaces, which is the null space of every non-negative combination of
them and so a property of the family rather than of a point.

## Usage

``` r
sum_struct_null_basis(components)
```

## Arguments

- components:

  A list of symmetric positive semidefinite matrices.

## Value

A matrix whose columns are an orthonormal basis, with zero columns when
the sum has full rank.

## Details

The components are stacked after each is divided by its own largest
entry. Without that normalization a component whose scale is many orders
below another's sinks below the tolerance and is read as absent, which
is the failure a rank taken from an assembled matrix already shows.
