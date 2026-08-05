# A Jet for One Variable

The jet of a scalar map of a single free value: its value and its four
derivatives sit in the pure components of that variable, everything else
being zero.

## Usage

``` r
jet_var(k, dv, lay)
```

## Arguments

- k:

  The index of the free value.

- dv:

  A list of five numbers: the value and four derivatives.

- lay:

  A layout from
  [`jet_layout`](https://statmodels7.github.io/parameters7/reference/jet_layout.md);
  taken from the jet itself unless given.

## Value

A jet.
