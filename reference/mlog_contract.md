# One Derivative Component by the Daleckii-Krein Contraction

Contracts the rotated directions of one index tuple against the
divided-difference table of the matching order, sums over **every
ordering** of the directions, and rotates the result back by \\Q\\.

## Usage

``` r
mlog_contract(tb, dirs)
```

## Arguments

- tb:

  The tables of
  [`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md),
  carrying at least `length(dirs) + 1` points.

- dirs:

  The free-value indices of the tuple, possibly repeated. Its length is
  the derivative order.

## Value

A symmetric numeric matrix of the side `tb$q` has, with no dimnames.

## Details

The sum runs over every ordering, never over the distinct ones, so a
tuple with repeated indices is counted with its multiplicity and needs
no correction afterwards. That is where a hand-written version of this
goes wrong: summing the distinct orderings and forgetting the \\\prod_j
m_j!\\ factor gives a result that is too small by exactly that factor,
which looks plausible.

The cost is the reason this family's higher orders are expensive: at
order 4 there are 24 orderings per component, each an \\O(p^4)\\
contraction. See
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
for the measurement.

## See also

[`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md)
for the input,
[`combinat_perms()`](https://statmodels7.github.io/parameters7/reference/combinat_perms.md)
for the orderings, and
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
for the representation.
