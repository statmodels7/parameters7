# The Subset Structure of One Derivative Order, Memoized

The index tuples of an order, the codes of every sub-multiset they
carry, and the codes of the sub-multisets of those, over which the inner
Leibniz rule sums. None of it depends on the free vector.

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

Memoized in an environment held on the object, because it depends on the
order and the number of free values alone. The cache is pure (the same
order always gives the same structure), so sharing it across the copies
that S7 makes of the object is safe.

## See also

[`ar_inv_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_inv_derivative.md),
the only caller.
