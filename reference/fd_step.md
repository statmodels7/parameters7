# Finite-Difference Step for a Free Value

The step for a central difference in one component of \\\eta\\, scaled
by the size of that component.

## Usage

``` r
fd_step(eta_k, order = 1L)
```

## Arguments

- eta_k:

  The value of the component.

- order:

  The derivative order the step is for.

## Value

A single positive number.

## Details

The free scale is unbounded, so unlike the response and parameter steps
of distributions7 this one has no boundary to be clamped away from: the
whole point of the unconstrained scale is that there is nowhere to fall
off, and `bounds` is left at its default for that reason.

The rule is numericals7's and is read from it rather than written out
again, so the step a fallback takes and the step the stencil library
documents cannot drift apart.

## See also

[`fd_step`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
