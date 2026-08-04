# Derivatives of the Logarithm of a Link

The derivatives of \\\log h(\eta)\\ up to fourth order, from the inverse
link's own derivatives.

## Usage

``` r
diag_dlog(link, eta, order)
```

## Arguments

- link:

  A linkfunctions7 link.

- eta:

  A numeric vector of free values.

- order:

  The derivative order, 1 to 4.

## Value

A numeric vector the length of `eta`.

## Details

Faa di Bruno's formula for the logarithm, written out: with \\u_m =
h^{(m)}/h\\ the successive orders are \\u_1\\, \\u_2 - u_1^2\\, \\u_3 -
3u_1u_2 + 2u_1^3\\ and \\u_4 - 4u_1u_3 - 3u_2^2 + 12u_1^2u_2 - 6u_1^4\\.
Dividing by \\h\\ once at the start keeps every term of order one in the
ratio rather than in the derivative, which matters at the ends of a
link's range.

## See also

[`diagonal_matrix`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`scaled_dlog`](https://statmodels7.github.io/parameters7/reference/scaled_dlog.md)
