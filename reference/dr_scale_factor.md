# The Scale Factor of a Derivative Component

Evaluates \\\partial^{S_D}(d_i d_j)\\ for every pair \\(i, j)\\, given
the multiset \\S_D\\ of scale indices. It is zero wherever \\S_D\\
contains an index other than \\i\\ or \\j\\, so the result is supported
on the rows and columns those indices name.

## Usage

``` r
dr_scale_factor(sd, tuple)
```

## Arguments

- sd:

  A 5 by `p` matrix of inverse-link derivatives, as returned by
  [`dr_scale_derivs`](https://statmodels7.github.io/parameters7/reference/dr_scale_derivs.md).

- tuple:

  An integer vector of scale indices, possibly empty and possibly with
  repeats.

## Value

A `p` by `p` numeric matrix.
