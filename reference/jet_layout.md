# The Bookkeeping a Jet Needs

The index tuples up to fourth order over \\d\\ variables, together with
the lookup from a tuple to its position, computed once and shared by
every jet in a calculation.

## Usage

``` r
jet_layout(d)
```

## Arguments

- d:

  The number of variables.

## Value

A list with `d`, the tuple lists `tuples`, and the environment `pos`
mapping a sorted tuple to its order and position.

## Details

The tuples are the package's own enumeration, so a jet's components are
keyed exactly as
[`param_tuple_names`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
keys a derivative list and nothing has to be reordered on the way out.

## See also

[`jet_var`](https://statmodels7.github.io/parameters7/reference/jet_var.md),
[`Ops.jet`](https://statmodels7.github.io/parameters7/reference/Ops.jet.md)

## Examples

``` r
lay <- jet_layout(2)
lay$d
#> [1] 2
length(lay$tuples[[2]])
#> [1] 3
```
