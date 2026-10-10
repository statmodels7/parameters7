# Finite-Difference Step for a Free Value

Returns the step for a central difference in one component of \\\eta\\,
scaled by the size of that component. It forwards to
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
at accuracy 2, so the step of a fallback here and the step that the
stencil library documents cannot drift apart.

## Usage

``` r
fd_step(eta_k, order = 1L)
```

## Arguments

- eta_k:

  The value of the component, a single number.

- order:

  The derivative order the step is for: 1, 2, 3 or 4.

## Value

A single positive number.

## Details

The rule is \\\varepsilon^{1/(k+2)} \max(1, \|\eta_k\|)\\ for a \\k\\-th
derivative, which balances the truncation error against the rounding
that the division by \\h^k\\ amplifies. In doubles that is \\6.1 \times
10^{-6}\\ at first order and \\1.2 \times 10^{-4}\\ at second, both
scaling up once \\\|\eta_k\|\\ passes 1.

The step is not clamped. The unconstrained scale has no boundary, so the
`bounds` argument of
[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
is left at its default.

## See also

[`numericals7::fd_step()`](https://statmodels7.github.io/numericals7/reference/fd_step.html)
for the rule, and
[`fd_along()`](https://statmodels7.github.io/parameters7/reference/fd_along.md),
which applies a stencil at this step.
