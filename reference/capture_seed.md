# Capture and Restore the Caller's Random Stream

`capture_seed()` reads `.Random.seed` from the global environment,
returning `NULL` when there is none, and `restore_seed()` puts back what
it was given. Together they let a function draw from a fixed seed, so
that its report is the same on two runs, without leaving the caller's
stream changed.

## Usage

``` r
capture_seed()

restore_seed(old)
```

## Arguments

- old:

  The value `capture_seed()` returned, or `NULL`.

## Value

`capture_seed()` returns the saved `.Random.seed`, an integer vector, or
`NULL`. `restore_seed()` returns `NULL` invisibly and is called for its
effect on the global environment.

## Details

A validator should report the same numbers on every run and should leave
the caller's random stream unchanged, and the two requirements conflict.
Drawing from the caller's stream makes the worst error reported vary
between runs; calling [`set.seed()`](https://rdrr.io/r/base/Random.html)
removes that variation but replaces the state of the caller's stream, so
a call inside a simulation changes the simulation. Saving the state on
entry and restoring it with
[`base::on.exit()`](https://rdrr.io/r/base/on.exit.html) gives the fixed
draw and leaves the caller's stream as it was.

`restore_seed(NULL)` removes `.Random.seed` again, which is the right
answer when the caller had never drawn a random number, so the seed that
the validator set is not left behind.

## Examples

``` r
# the property the pair exists for, through the public interface: a
# validator call leaves the caller's stream exactly as it found it
set.seed(42)
want <- rnorm(1)
set.seed(42)
invisible(check_parameter(log_cholesky(2), verbose = FALSE))
identical(rnorm(1), want)
#> [1] TRUE
```
