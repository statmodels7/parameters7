# Divided Differences of the Exponential

Returns \\e\[\lambda_1, \dots, \lambda_m\]\\ by the Opitz theorem: build
the upper bidiagonal matrix with the arguments on the diagonal and ones
above it, exponentiate it, and read the divided difference off the
top-right corner. The value is exact under repeated arguments, where the
recursive quotient \\(f\[\lambda_2..\] - f\[\lambda_1..\])/(\lambda_m -
\lambda_1)\\ divides by zero, and accurate under nearly repeated ones,
where the quotient cancels.

## Usage

``` r
dd_exp(lams)
```

## Arguments

- lams:

  The arguments, which may be in any order and may repeat. A divided
  difference is symmetric in its arguments, so the order does not
  matter. The number of arguments is at most 5, which is the number that
  fourth-order derivatives need.

## Value

A single number.

## Details

When two arguments nearly coincide, the quotient loses accuracy as the
gap shrinks, while the Opitz value stays accurate to rounding.
Near-repeated eigenvalues are ordinary here: any \\S\\ with a symmetry
has them, and a scalar multiple of the identity has an exactly repeated
pair.

## See also

[`mlog_expm_small()`](https://statmodels7.github.io/parameters7/reference/mlog_expm_small.md),
which does the exponential, and
[`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md),
which builds whole tables of these.
