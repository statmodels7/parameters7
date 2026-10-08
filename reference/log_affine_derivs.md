# Derivatives of a Sum of Logarithms of Affine Functions

Returns the derivatives of orders 1 to `order` at \\r\\ of \\\sum_t c_t
\log(a_t + b_t r)\\, using

\$\$\frac{\mathrm{d}^k}{\mathrm{d}r^k}\log(a + br) =
(-1)^{k-1}(k-1)!\\\frac{b^k}{(a + br)^k}.\$\$

## Usage

``` r
log_affine_derivs(r, terms, order)
```

## Arguments

- r:

  The point, a single number.

- terms:

  A list of numeric triples `c(coefficient, a, b)`, one per logarithm.

- order:

  The highest order wanted, an integer from 1 to 4.

## Value

A list of `order` numbers, the first to the `order`-th derivative.

## Details

Both two-value families have a log-determinant of this shape: compound
symmetry from its two distinct eigenvalues, \\\log\\1+(p-1)\rho\\ +
(p-1)\log(1-\rho)\\, and AR(1) from \\\|R\| = (1-\rho^2)^{p-1}\\, which
factors as \\(p-1)\\\log(1+\rho) + \log(1-\rho)\\\\. Writing the
derivatives once here keeps them from being transcribed twice, which is
where a sign would be lost.

The result is in the **correlation**, not in the free value; the caller
chains it onto the link with
[`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md).

## See also

[`cs_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_terms.md)
and
[`ar1_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_terms.md)
for the two term lists, and
[`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md),
which chains the result onto the link.
