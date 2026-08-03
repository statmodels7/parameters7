# Free Vectors to Sweep a Structure Over

A set of free vectors covering the ordinary range and a spread one.

## Usage

``` r
sweep_etas(s, n = 4L)
```

## Arguments

- s:

  A
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  object.

- n:

  The number of random vectors.

## Value

A list of numeric vectors.

## Details

The spread of the last vector is deliberately moderate. The free scale
is unbounded, so no value of \\\eta\\ is inadmissible, but the matrix
built from widely separated free values can be singular in double
precision – spreading a log-Cholesky structure over twenty-eight units
of log gives a condition number around \\10^{28}\\ – and a comparison
that fails there is a statement about the arithmetic rather than about
the structure. The scaling that a rank-deficient family must survive is
a property of its components rather than of a point, so it is tested
where it arises, against the declared null space, and not by driving
every family off the edge of double precision.
