# Derivatives of the Logarithm of a Link

Returns \\\mathrm{d}^m \log h(\eta) / \mathrm{d}\eta^m\\ for \\m\\ up to
four, assembled from the inverse link's own derivatives. Every
log-determinant derivative of a diagonal or scaled family is one of
these, multiplied by a count of entries.

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

  The derivative order: 1, 2, 3 or 4.

## Value

A numeric vector the length of `eta`.

## Details

Faa di Bruno's formula for the logarithm, written out. With \\u_m =
h^{(m)}/h\\ the four orders are

\$\$u_1, \quad u_2 - u_1^2, \quad u_3 - 3u_1u_2 + 2u_1^3, \quad u_4 -
4u_1u_3 - 3u_2^2 + 12u_1^2u_2 - 6u_1^4.\$\$

The division by \\h\\ happens once, at the start, so every term is a
ratio of order one, never a derivative that can be large on its own.
That matters at the ends of a link's range, where \\h\\ and \\h^{(m)}\\
can both be extreme while the ratio is ordinary.

Under the log link \\h = e^\eta\\ gives \\u_m = 1\\ for every \\m\\, so
the first order is 1 and the second, third and fourth are exactly 0,
which is why a default
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
reports a linear log-determinant.

## See also

[`diag_logdet_higher()`](https://statmodels7.github.io/parameters7/reference/diag_logdet_higher.md),
[`param_d3logdet.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.ScaledMatrixParam.md)
and
[`param_d4logdet.ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.ScaledMatrixParam.md),
the callers, and
[`linkfunctions7::dlinkinv()`](https://statmodels7.github.io/linkfunctions7/reference/dlinkinv.html)
for the link derivatives read.
