# The Blocks' Derivative Components, Keyed by Local Index Tuple

Fetches one block's derivatives of a given order and re-keys them by the
**sorted local index tuple**, as `"1"`, `"1,2"`, `"2,2"` and so on.

## Usage

``` r
block_derivs_by_tuple(block, eta, order)
```

## Arguments

- block:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  one block of the composite.

- eta:

  The block's own stretch of the free vector, of length `block@n_free`.

- order:

  The derivative order: 1, 2, 3 or 4.

## Value

A list of the block's own derivative matrices of that order, named by
the sorted local index tuples.

## Details

The re-keying lets a lookup from the composite's enumeration succeed
without assuming the two enumerations correspond. A block's own
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
are built from its own `free_names`, which the composite has prefixed
with a label, so the two name sets differ; the sorted index tuple is the
one key both sides can compute. Sorting matters because a derivative is
symmetric in its indices and the two enumerations need not order a tuple
the same way.

## See also

[`block_diag_derivs()`](https://statmodels7.github.io/parameters7/reference/block_diag_derivs.md),
the only caller.
