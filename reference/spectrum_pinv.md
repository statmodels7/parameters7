# Moore-Penrose Inverse From a Parameter's Spectrum

Builds the pseudo-inverse of \\M(\eta)\\ as \\\sum\_{j \in
\mathrm{keep}} \lambda_j^{-1} v_j v_j^\top\\, over the directions
[`param_spectrum()`](https://statmodels7.github.io/parameters7/reference/param_spectrum.md)'s
`keep` flag marks. It is what the rank-deficient branches of the
base-class log-determinant derivatives use in place of \\M^{-1}\\ in the
trace identities.

## Usage

``` r
spectrum_pinv(sp)
```

## Arguments

- sp:

  The result of
  [`param_spectrum()`](https://statmodels7.github.io/parameters7/reference/param_spectrum.md):
  a list with `values`, `vectors` and `keep`.

## Value

A symmetric numeric matrix of the same side as `sp$vectors`. When `keep`
is all `FALSE`, which happens only at rank 0, the zero matrix.

## See also

[`param_spectrum()`](https://statmodels7.github.io/parameters7/reference/param_spectrum.md)
for the input, and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md),
which rejects a deficient family rather than handing this back.
