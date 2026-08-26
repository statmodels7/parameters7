# Validate a Free Vector Against a Parameter

Checks that `eta` is a finite numeric vector of the length the parameter
declares, and returns it with its names stripped. Called in the body of
every generic before dispatch, so a parameter written outside the
package inherits the check without doing anything.

## Usage

``` r
check_eta(s, eta)
```

## Arguments

- s:

  A
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  object, whose `n_free` and `free_names` are read.

- eta:

  The free vector supplied by the caller.

## Value

`eta`, as an unnamed numeric vector of length `s@n_free`.

## Details

Three conditions are checked, each with its own message. A non-numeric
`eta` throws `'eta' must be numeric.`; a wrong length throws a message
naming both counts and listing the family's `free_names`, so a caller
who has mismatched two parametrizations can see which is which; and any
`NA`, `NaN` or infinite entry throws
`'eta' must be finite: the free scale has no boundary to reach.` The
last message states the reason a non-finite free value is a caller
error: the unconstrained scale has no edge, so an infinity there is a
runaway rather than a limit reached.

Names are stripped for the reason `align_theta()` strips them in
distributions7. A value that has been through a link comes back carrying
its own name, which means nothing on a number and would otherwise appear
in the dimnames of the matrix built from it.
