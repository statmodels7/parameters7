# Default Log-Determinant Hessian

The method every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
inherits when it registers no
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
of its own. It evaluates the identity that follows from differentiating
[`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)'s
trace once more, using \\\partial_l M^{-1} = -M^{-1}(\partial_l
M)M^{-1}\\,

\$\$\partial\_{kl} \log\|M\| = \mathrm{tr}\\\left(M^{+} \partial\_{kl}
M\right) - \mathrm{tr}\\\left(M^{+} (\partial_k M)\\ M^{+} (\partial_l
M)\right),\$\$

with \\M^{+}\\ the pseudo-inverse over the directions the declared rank
keeps. The second trace is formed as `sum(t(mi %*% dk) * (mi %*% dl))`,
which is the trace of the product without the product being multiplied
out.

Exact given the derivative arrays, so the accuracy is theirs. With
analytic arrays the answer agrees with a closed form to \\1 \times
10^{-14}\\; with numerical ones, to \\6 \times 10^{-8}\\.

## Arguments

- s:

  A
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A numeric vector of `choose(s@n_free + 1, 2)` entries, keyed as
`param_tuple_names(s)` and in that order.

## See also

[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
for the generic,
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
for the arrays this reads, and
[`param_d3logdet.matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.matrix_parameter.md)
for the order above, which differences this one.
