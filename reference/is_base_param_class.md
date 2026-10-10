# Test for the Package's Own Base Classes

Returns whether an S7 class is
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
or
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
the two base classes on which this package registers its fallback
methods. The test distinguishes a method that a subclass wrote from one
that it inherited, which is what
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reports and what
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
acts on.

## Usage

``` r
is_base_param_class(cls)
```

## Arguments

- cls:

  An S7 class object, as `S7::S7_class(x)` returns or as
  `attr(method, "signature")[[1]]` holds. Anything else returns `FALSE`
  without an error, [`attr()`](https://rdrr.io/r/base/attr.html) on a
  non-class giving `NULL`.

## Value

`TRUE` when `cls` is
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
or
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
`FALSE` for any subclass and for anything that is not an S7 class.

## Details

Two comparisons, in order. Identity is tried first, since it is the
usual case and costs nothing. Then the class name together with the
owning package, because identity does not survive a package's code being
re-evaluated instead of loaded: a class rebuilt from the same definition
is a different object. Coverage tools re-evaluate, so an identity-only
test passes every ordinary check and fails under `covr` alone.

## See also

[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md),
which calls this on the class a method was registered against.
