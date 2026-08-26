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

A validator wants both properties and they pull against each other.
Drawing from the caller's stream makes the worst error reported move
between runs; calling [`set.seed()`](https://rdrr.io/r/base/Random.html)
fixes that and replaces whatever state the caller had, so a call inside
a simulation silently changes the simulation. Saving the state on entry
and restoring it with
[`base::on.exit()`](https://rdrr.io/r/base/on.exit.html) gives the fixed
draw and leaves the caller alone.

`restore_seed(NULL)` removes `.Random.seed` again, which is the right
answer when the caller had never drawn a random number: leaving the seed
the validator set behind would be the leak this exists to prevent.

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
