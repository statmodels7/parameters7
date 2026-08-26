# The Eigendecomposition and Divided-Difference Tables at a Point

Everything the derivative contractions need: the eigendecomposition of
\\S\\, the rotated basis directions, and the divided-difference tables
of the orders asked for, computed once per free vector.

## Usage

``` r
mlog_tables(s, eta, order)
```

## Arguments

- s:

  A
  [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`.

- order:

  The highest derivative order wanted: 1, 2, 3 or 4.

## Value

A list with `q` and `lam`, the eigenvectors and eigenvalues of \\S\\;
`e`, a list of `s@n_free` rotated directions; and the divided-difference
tables `dd2` up to `dd<order+1>`.

## Details

Everything here depends on \\\eta\\ alone, never on the index tuple, so
it is computed once and every component of the order reads it. The
rotated directions are \\Q^\top E_k Q\\ for each basis direction
\\E_k\\, and the tables hold one divided difference per combination of
eigenvalues with repetition, up to `order + 1` points.

## See also

[`dd_exp()`](https://statmodels7.github.io/parameters7/reference/dd_exp.md)
for the entries of the tables, and
[`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md),
which reads them.
