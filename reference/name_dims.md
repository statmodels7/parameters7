# Name the Rows and Columns of a Parameter's Matrix

Applies the dimension labels every matrix a parameter produces carries:
`"v1"`, `"v2"`, ..., `"vp"` on both margins, `p` being the parameter's
`dimension`. One convention across the families, so a consumer can read
a printed covariance without knowing which parametrization built it.

It covers the value and the four derivative orders. A factor and a solve
are left as each family produces them, the families disagreeing there,
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
labels its own margins `"s1"`, `"s2"`, ..., its rows being states rather
than variables.
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)'s
shapes-and-names check reads the margins of
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
so the convention is tested rather than only stated.

## Usage

``` r
name_dims(m, s)
```

## Arguments

- m:

  A numeric matrix, `p` by `p`. Not checked; the callers are the
  package's own
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  methods and the helpers assembling the derivative arrays.

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, whose `dimension` supplies `p`.

## Value

`m`, with `dimnames` set to `list(v1..vp, v1..vp)`. Any dimnames it
arrived with are replaced.
