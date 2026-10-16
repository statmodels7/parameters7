# Test for the Package's Own Base Classes

Returns whether an S7 class is
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
or
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
the two base classes on which this package registers its fallback
methods. The test distinguishes a method that a subclass wrote from one
that it inherited, and
[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
reports this distinction while
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
acts on it.

## Usage

``` r
is_base_param_class(cls)
```

## Arguments

- cls:

  An S7 class object, as `S7::S7_class(x)` returns or as
  `attr(method, "signature")[[1]]` holds. Anything else returns `FALSE`
  without an error, because [`attr()`](https://rdrr.io/r/base/attr.html)
  on a non-class gives `NULL`.

## Value

`TRUE` when `cls` is
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
or
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
`FALSE` for any subclass and for anything that is not an S7 class.

## Details

The test applies two comparisons in order. It tries identity first,
since this is the usual case and costs nothing. It then compares the
class name together with the owning package, because identity does not
survive when a package's code is re-evaluated instead of loaded: a class
rebuilt from the same definition is a different object. Coverage tools
re-evaluate, so an identity-only test passes every ordinary check and
fails only under `covr`.

## See also

[`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md),
which calls this on the class against which a method was registered.
