# The Ordered-Partition Sum for the Derivatives of an Inverse

Computes the derivatives of one order of \\N = M^{-1}\\ from the
derivatives of \\M\\, as the sum over the ordered set partitions of the
index tuple of \\(-1)^q N A\_{B_1} N \cdots A\_{B_q} N\\. It is the
default route of
[`param_inv_d3()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
and
[`param_inv_d4()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md),
and the route of
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
over a family that has no exact method.

## Usage

``` r
inv_partition_sum(inner, eta, order)
```

## Arguments

- inner:

  An object inheriting from
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  of full rank.

- eta:

  A numeric vector of length `inner@n_free`, already checked.

- order:

  The derivative order, 1 to 4.

## Value

An unnamed list of `inner@dimension` square matrices without dimnames,
ordered as `param_tuple_indices(inner, order)`.

## Details

The derivative arrays of \\M\\ of every order up to `order` are fetched
once and keyed by their sorted index tuple, so that a block appearing in
many partitions is read and not recomputed.

## See also

[`ordered_set_partitions()`](https://statmodels7.github.io/parameters7/reference/ordered_set_partitions.md)
for the enumeration and
[`chol_inv_dk()`](https://statmodels7.github.io/parameters7/reference/chol_inv_dk.md)
for the exact route of a log-Cholesky parameter.
