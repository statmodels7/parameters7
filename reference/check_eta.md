# Validate a Free Vector Against a Parameter

Checks that `eta` is a finite numeric vector of the length the parameter
declares, and returns it unnamed.

## Usage

``` r
check_eta(s, eta)
```

## Arguments

- s:

  A
  [`parameter`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object.

- eta:

  The free vector supplied by the caller.

## Value

`eta`, as an unnamed numeric vector.

## Details

Called in the body of every generic before dispatch, so that a parameter
written outside the package inherits the check without doing anything.
The names are stripped for the reason `align_theta()` strips them in
distributions7: a value that has been through a link comes back carrying
its own name, which is meaningless on a number and would leak into the
dimnames of the result.
