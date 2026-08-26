# The Column of a Packed Derivative Record

Locates a derivative component in a row of
[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)'s
output. The value sits in column 1 and the full tensors of orders one to
four follow it in row-major order, so order \\k\\ begins at column \\1 +
\sum\_{j\<k} n^{j}\\ and the tuple \\(t_1, \dots, t_k)\\ sits \\\sum_i
(t_i - 1) n^{k-i}\\ further on.

## Usage

``` r
ar_pack_col(n, order = 0L, tuple = NULL)
```

## Arguments

- n:

  The number of free values, \\q + 1\\.

- order:

  The derivative order 1 to 4, or 0 for the value.

- tuple:

  The index tuple, 1-based, of length `order`. Ignored at order 0.

## Value

A single column index.

## Details

At \\n = 3\\ the value is column 1, the three first-order components are
columns 2 to 4, and the nine second-order ones are columns 5 to 13:
\\(1,1)\\ is 5, \\(1,2)\\ is 6 and \\(3,3)\\ is 13. The tensor is stored
**full**, with one column per ordered tuple, so \\(1,2)\\ and \\(2,1)\\
are two columns holding the same number; the caller reads whichever one
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
hands it.

## See also

[`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
for the layout this describes.
