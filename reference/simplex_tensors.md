# Derivative Tensors of the Softmax Map

Computes the arrays \\D_1\[a, b\]\\ to \\D_4\[a, b, c, d, e\]\\ of
derivatives of \\\pi\\ in the free values, by applying the product rule
to \\\partial_b \pi_a = \pi_a(\delta\_{ab} - \pi_b)\\ once for each
order up to `order`. The first index runs over the \\K\\ categories of
the value, the rest over the \\K - 1\\ free values.

## Usage

``` r
simplex_tensors(pi_full, order = 4L)
```

## Arguments

- pi_full:

  The value, a vector of length \\K\\ summing to 1, as
  [`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
  returns it.

- order:

  The highest order to compute: 1, 2, 3 or 4.

## Value

A list with elements `d1` to `d<order>`: `d1` is `K` by `K-1`, `d2` is
`K` by `K-1` by `K-1`, and so on. In each array the first index runs
over the categories of the value and the rest run over the free values.

## Details

Each array is one more application of the product rule to the one below,
so the cost of the fourth is \\K(K-1)^4\\ entries and the whole set is
built in one pass. The entries follow from the product rule, and the
construction does not use finite differences.

Every slice sums to zero over its first index, because \\\sum_a \pi_a\\
is the constant 1. This is the identity that
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
tests.

## See also

[`simplex_components()`](https://statmodels7.github.io/parameters7/reference/simplex_components.md),
which slices these into the named lists that the generics return, and
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
for the value from which they are built.
