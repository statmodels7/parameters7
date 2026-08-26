# The Softmax Point Behind a Free Vector

Computes the value \\\pi\\ from a free vector, appending the reference
category's implicit 0 and taking the softmax through the **log-sum-exp
shift**: the largest entry is subtracted before exponentiating, so a
large free value saturates instead of overflowing.

## Usage

``` r
simplex_point(eta)
```

## Arguments

- eta:

  A numeric vector of length \\K - 1\\. Not checked; the callers have
  been through the generic.

## Value

A numeric vector of length \\K\\ with non-negative entries summing to 1,
and no names. The names `p1` ... `pK` are applied by
[`param_value.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_value.SimplexParam.md).

## Details

The naive expression \\e^{\eta_a}/(1 + \sum_b e^{\eta_b})\\ is `Inf/Inf`
from about \\\eta = 710\\. With the shift, a free value of 800 returns
\\\pi_1 = 1\\ and \\\pi_K = 0\\ with the vector summing to exactly 1,
and at 500 it returns \\\pi_K = 7 \times 10^{-218}\\, still
representable. The value is then on the boundary of the simplex, which
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
cannot invert.

## See also

[`param_value.SimplexParam()`](https://statmodels7.github.io/parameters7/reference/param_value.SimplexParam.md)
and
[`param_value.TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/param_value.TransitionMatrixParam.md),
the two callers, and
[`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
for the derivatives of the same map.
