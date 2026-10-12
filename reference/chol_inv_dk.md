# Higher Derivatives of a Log-Cholesky Inverse

Computes the derivatives of one order of \\M^{-1} = G^\top G\\, with \\G
= L^{-1}\\, by the Leibniz rule on \\L G = I\\ and on \\G^\top G\\, as
[`param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
describes. The derivatives of \\G\\ are kept by sorted multiset, so that
a multiset shared by several components is computed once.

## Usage

``` r
chol_inv_dk(s, eta, order)
```

## Arguments

- s:

  A
  [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  object.

- eta:

  Its free vector, already checked.

- order:

  The derivative order, 3 or 4 in the package's use; any positive order
  is computed the same way.

## Value

A list with one matrix per tuple in `param_tuple_indices(s, order)`,
named as `param_tuple_names(s, order)`.

## See also

[`chol_inv_pieces()`](https://statmodels7.github.io/parameters7/reference/chol_inv_pieces.md)
for the first two orders and
[`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md)
for the Leibniz sum of a Gram product.
