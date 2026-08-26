# Divided Differences of the Exponential

Returns \\e\[\lambda_1, \dots, \lambda_m\]\\ by the **Opitz theorem**:
build the upper bidiagonal matrix with the arguments on the diagonal and
ones above it, exponentiate it, and read the divided difference off the
top-right corner. Exact under repeated arguments, where the recursive
quotient \\(f\[\lambda_2..\] - f\[\lambda_1..\])/(\lambda_m -
\lambda_1)\\ divides by zero, and accurate under nearly repeated ones,
where it cancels.

## Usage

``` r
dd_exp(lams)
```

## Arguments

- lams:

  The arguments, in any order and with repeats allowed. A divided
  difference is symmetric in its arguments, so the order does not
  matter. At most 5 of them, the number fourth-order derivatives need.

## Value

A single number.

## Details

Measured at \\\lambda = (0.5,\\ 0.5 + g)\\, against the exact limit
\\e^{0.5} = 1.648721270700127\\ reached at \\g = 0\\:

|             |                   |                   |
|-------------|-------------------|-------------------|
| gap \\g\\   | Opitz             | quotient          |
| \\10^{-3}\\ | 1.649545906191066 | 1.649545906190929 |
| \\10^{-6}\\ | 1.648722095061039 | 1.648722095023754 |
| \\10^{-9}\\ | 1.648721271524488 | 1.648721206226264 |

The quotient has lost eight digits at a gap of \\10^{-9}\\ and the Opitz
value one. Near-repeated eigenvalues are ordinary here: any \\S\\ with a
symmetry has them, and an exactly repeated pair is what a scalar
multiple of the identity gives.

## See also

[`mlog_expm_small()`](https://statmodels7.github.io/parameters7/reference/mlog_expm_small.md),
which does the exponential, and
[`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md),
which builds whole tables of these.
