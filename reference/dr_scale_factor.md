# The Scale Factor of a Derivative Component

Evaluates \\\partial^{S_D}(d_i d_j)\\ for every pair \\(i, j)\\ at once,
given the multiset \\S_D\\ of scale indices. It is zero wherever \\S_D\\
contains an index naming neither \\i\\ nor \\j\\, so the result is
**supported on the rows and columns those indices name**, and that
sparsity is half of why the composition costs so little.

## Usage

``` r
dr_scale_factor(sd, tuple)
```

## Arguments

- sd:

  A 5 by \\p\\ matrix of inverse-link derivatives, as returned by
  [`dr_scale_derivs()`](https://statmodels7.github.io/parameters7/reference/dr_scale_derivs.md).

- tuple:

  An integer vector of scale indices, possibly empty and possibly with
  repeats. Its length is the number of scale indices in the component,
  which is at most the derivative order.

## Value

A symmetric \\p\\ by \\p\\ numeric matrix.

## Details

Three cases, and the arithmetic differs in each:

- **no indices.** The factor is \\d_i d_j\\, an outer product, and that
  is the case a component differentiating the correlation alone falls
  into.

- **one distinct index \\k\\, with multiplicity \\m\\.** Off the
  diagonal the entry carries one factor of \\d_k\\, so it is \\d_k^{(m)}
  d_j\\; on the diagonal it carries two, so it is the Leibniz expansion
  \\\sum_r \binom{m}{r} d_k^{(r)} d_k^{(m-r)}\\ of \\\partial^m d_k^2\\.

- **two distinct indices \\i\\ and \\j\\.** Only the entries \\(i,j)\\
  and \\(j,i)\\ survive, each \\d_i^{(m_i)} d_j^{(m_j)}\\.

Three or more distinct indices give the zero matrix: an entry of
\\\Sigma\\ carries at most two scales, so a third differentiation in a
new scale annihilates it.

## See also

[`dr_prod_derivs()`](https://statmodels7.github.io/parameters7/reference/dr_prod_derivs.md),
the only caller, and
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
for the factorization this is half of.
