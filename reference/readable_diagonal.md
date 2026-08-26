# Quantities That Are Separate Links of Separate Free Values

Assembles a
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
declaration for a family whose quantities are one scalar link each of
one free value each. The Jacobian is then **diagonal**, its \\k\\-th
entry the inverse link's first derivative at \\\eta_k\\, and there is
nothing else to compute.

[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
and
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
all declare through this, which is why their three methods are one call
each.

## Usage

``` r
readable_diagonal(links, eta, nm, transform, label)
```

## Arguments

- links:

  A list of linkfunctions7 links, one per quantity, in the order of the
  free values they read. `links[[k]]` must be the link of `eta[k]`; the
  correspondence is positional and is not checked.

- eta:

  A numeric vector of free values, at least as long as `links`.

- nm:

  A character vector naming the quantities, the same length as `links`.

- transform:

  A character vector naming the scale each interval is built on, the
  same length as `nm`.

- label:

  A single string naming the block.

## Value

A list with `value`, `jacobian`, `transform` and `label`, as
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
describes. The Jacobian is `length(nm)` by `length(eta)` with `nm` as
its row names, and is diagonal.

## See also

[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
for the contract, and
[`param_readable.AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/param_readable.AutoregressiveParam.md)
for the one family whose Jacobian is not diagonal.
