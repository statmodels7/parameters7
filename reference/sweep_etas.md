# Free Vectors to Sweep a Parameter Over

Builds the set of free vectors
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
runs its battery at: the origin, two constant vectors at \\\pm 0.5\\,
`n` drawn uniformly from \\(-1.5, 1.5)\\, and one evenly spread from
\\-2\\ to \\2\\. Every check reports the worst discrepancy over the
whole set, so a family that is right at the origin and wrong away from
it is caught.

## Usage

``` r
sweep_etas(s, n = 4L)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, whose `n_free` sets the length.

- n:

  The number of random vectors, defaulting to 4. The set returned has
  `n + 4` members, or one member, `numeric(0)`, for a family with
  nothing to estimate.

## Value

A list of numeric vectors, each of length `s@n_free`.

## Details

The spread of the last vector is deliberately moderate. The free scale
is unbounded, so no \\\eta\\ is inadmissible, but the matrix built from
widely separated free values can be singular in double precision:
spreading a log-Cholesky parameter over twenty-eight units of log gives
a condition number around \\10^{28}\\, and a comparison that fails there
says something about the arithmetic and nothing about the
parametrization.

The scaling a rank-deficient family has to survive is a property of its
components, the same at every point, so it is tested where it arises,
against the declared null space, instead of by driving every family off
the edge of double precision.

The `n` random vectors are drawn from whatever random stream is current.
No seed is set here:
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
the only caller, fixes one before calling and restores the caller's own
on exit, so the report is reproducible and the caller's stream is left
as it was found.

## See also

[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
the only caller.
