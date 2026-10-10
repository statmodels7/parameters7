# A Link Inverse and Its Derivatives, to a Given Order

Returns \\g^{-1}(\eta)\\ and its derivatives of orders 1 to `order`,
from which the families built on a link compose their maps. Each order
is one call of the link's own generic (`dlinkinv()` to `d4linkinv()`),
and the orders above `order` are not evaluated.

## Usage

``` r
linkinv_upto(link, e, order)
```

## Arguments

- link:

  A linkfunctions7 link.

- e:

  A numeric vector of free values.

- order:

  The highest order wanted, an integer from 0 to 4.

## Value

A list of `order + 1` elements: the value, then the derivatives in
order, each the shape of `e`.

## See also

[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md),
which composes onto these.
