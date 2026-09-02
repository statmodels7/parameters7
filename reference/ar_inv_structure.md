# The Subset Structure of One Derivative Order, Memoized

The index tuples of an order, the codes of every sub-multiset they
carry, and the codes of the sub-multisets of THOSE, which the inner
Leibniz sums over. None of it depends on the free vector.

## Usage

``` r
ar_inv_structure(s, order)
```

## Arguments

- s:

  An
  [`AutoregressiveInvParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveInvParam.md)
  object.

- order:

  The derivative order, 1 to 4.

## Value

A list with `cd` and `sub`, both as
[`ar_inv_codes()`](https://statmodels7.github.io/parameters7/reference/ar_inv_codes.md)
returns, and `names`, the component names of the order.

## Details

Memoized in an environment held on the object, because it is a function
of the order and the number of free values alone and was 59 per cent of
a fourth derivative when rebuilt at every call. The cache is pure – the
same order always gives the same structure – so sharing it across the
copies S7 makes of the object is safe.

## See also

[`ar_inv_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_inv_derivative.md),
the only caller.
