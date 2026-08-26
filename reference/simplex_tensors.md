# Derivative Tensors of the Softmax Map

The arrays \\D_1\[a, b\]\\ to \\D_4\[a, b, c, d, e\]\\ of derivatives of
\\\pi\\ in the free values, built by applying the product rule to
\\\partial_b \pi_a = \pi_a(\delta\_{ab} - \pi_b)\\ as many times as the
order asks. The first index runs over the \\K\\ categories of the value,
the rest over the \\K - 1\\ free values.

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

  The highest order wanted: 1, 2, 3 or 4.

## Value

A list with elements `d1` to `d<order>`: `d1` is `K` by `K-1`, `d2` is
`K` by `K-1` by `K-1`, and so on, the first index running over the
categories of the value and the rest over the free values.

## Details

Each array is one more application of the product rule to the one below,
so the cost of the fourth is \\K(K-1)^4\\ entries and the whole set is
built in one pass. Nothing is differenced.

Every slice sums to zero over its first index, \\\sum_a \pi_a\\ being
the constant 1, which is the identity
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
tests.

## See also

[`simplex_components()`](https://statmodels7.github.io/parameters7/reference/simplex_components.md),
which slices these into the named lists the generics return, and
[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
for the value they are built from.
