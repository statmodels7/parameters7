# Negate Without Changing the Container

Returns `-x` for a numeric vector and the elementwise negation for a
list, so that the log-determinant derivatives of
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
come back in the shape the inner family returned them in.

## Usage

``` r
neg_like(x)
```

## Arguments

- x:

  A numeric vector or a list of numbers.

## Value

`x` negated, of the same type.

## Details

The container is part of the contract:
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
subtracts the result of
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
from a reference, which a list does not support, so an
[`lapply()`](https://rdrr.io/r/base/lapply.html) over a numeric vector
turns a correct value into an error several frames away.

## See also

[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md),
the only caller.
