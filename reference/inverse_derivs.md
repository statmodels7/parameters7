# Derivatives of an Inverse Parameter of a Given Order

The ordered-block-partition sum of
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md),
run over the enumeration the class is keyed by.

## Usage

``` r
inverse_derivs(s, eta, order)
```

## Arguments

- s:

  An
  [`InverseParam()`](https://statmodels7.github.io/parameters7/reference/InverseParam.md)
  object.

- eta:

  A numeric vector of length `s@n_free`, already checked by the generic.

- order:

  The derivative order, 1 to 4.

## Value

A named list of `s@dimension` square matrices, keyed and ordered as
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
says.

## Details

The inner family's derivative arrays of every order up to `order` are
fetched once and keyed by their sorted index tuple, so a block appearing
in many partitions is read and not recomputed.

## See also

[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
for the formula.
