# Third and Fourth Derivatives of a Diagonal Matrix

Assembles a whole derivative order of a diagonal family, orders three
and four. A diagonal family is separable, so a component is the zero
matrix unless every index of the tuple names the **same** free value,
and a surviving one carries \\h'''(\eta_k)\\ or \\h''''(\eta_k)\\ in the
entries that value owns. Most of the list is therefore exactly zero: at
\\p = 3\\ only 3 of the 10 third-order components are non-zero.

## Usage

``` r
diag_higher(s, eta, order)
```

## Arguments

- s:

  A
  [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The derivative order: 3 or 4.

## Value

A list of `choose(s@n_free + order - 1, order)` diagonal matrices keyed
as `param_tuple_names(s, order)` and in that order, each `s@dimension`
by `s@dimension` with dimnames `v1`, `v2`, ...

## See also

[`diag_owner()`](https://statmodels7.github.io/parameters7/reference/diag_owner.md)
for the ownership map, and
[`param_d3.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3.DiagMatrixParam.md)
and
[`param_d4.DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.DiagMatrixParam.md),
the two callers.
