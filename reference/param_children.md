# The Parameters Held Inside a Parameter

Returns the
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
objects found in `s@param_params`, at any depth of nesting in lists: the
blocks of a
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
the inner parameter of a
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
or an
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md),
the terms of a
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md).

## Usage

``` r
param_children(s)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

## Value

A list of
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
objects, empty for a family that holds none.

## See also

[`param_analytic()`](https://statmodels7.github.io/parameters7/reference/param_analytic.md),
which reads them.
