# Construct a Simplex Parameter

Returns an object holding the additive log-ratio map onto the open
simplex: a probability vector of length \\K\\ from \\K - 1\\
unconstrained free values, with the last category as reference,

\$\$\pi_a = \frac{e^{\eta_a}}{1 + \sum_b e^{\eta_b}}, \qquad \pi_K =
\frac{1}{1 + \sum_b e^{\eta_b}}.\$\$

Every free vector gives positive entries summing to exactly 1, so a
mixture weight, a categorical probability or a latent-state distribution
can be estimated without a constraint. This is the first family here
whose value is **not a matrix**; see
[`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
for what that costs.

## Usage

``` r
simplex(n_cat)
```

## Arguments

- n_cat:

  The number of categories \\K\\, **at least 2**. A single integer; `1`,
  a fraction, `NA` and a vector all throw
  `'n_cat' must be a single integer of at least 2.` A one-category
  simplex is the constant 1 and has nothing to estimate.

## Value

An object of class
[`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md),
with `n_free` equal to `n_cat - 1`, `free_names` `alr1` ... `alr(K-1)`,
`param_name` `"simplex"`, and `param_params` holding `n_cat`. It carries
no `dimension`, `rank`, `null_basis` or `role`, being a
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
and never a
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md).

## The derivatives close over the value

With \\Z = 1 + \sum_b e^{\eta_b}\\, each \\\pi_a\\ for \\a \< K\\ is a
derivative of \\\log Z\\, so the derivative tensors of \\\pi\\ are the
cumulants of a categorical indicator and every order follows from one
rule applied repeatedly:

\$\$\partial_b \pi_a = \pi_a(\delta\_{ab} - \pi_b).\$\$

At first order that is the covariance matrix of a categorical indicator,
which is why the tensors are the same objects a multinomial score
already carries. All four orders are closed form and no stencil is used.

## An identity worth checking against

\\\sum_a \pi_a = 1\\ at every \\\eta\\, so differentiating it gives
\\\sum_a \partial \pi_a = 0\\, and the same at every higher order. Every
derivative component therefore sums to zero over the value index.
Measured over all four orders the worst sum is \\1.7 \times 10^{-17}\\,
and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
runs exactly this check: a derivative array that does not sum to zero is
wrong whatever else it agrees with.

## Large free values saturate instead of overflowing

[`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
applies the log-sum-exp shift, so a free value of 800 gives \\\pi_1 =
1\\ and \\\pi_K = 0\\ with the vector still summing to exactly 1, where
the naive expression would divide `Inf` by `Inf`. Note that the value
then sits on the **boundary** of the simplex, so
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
cannot invert it: \\\log(0)\\ is not finite.

## Why not stick-breaking

The other common chart chains through \\K - 1\\ nested logistic maps, so
its derivatives compose that many times and are not symmetric in the
categories. The additive log-ratio's close over \\\pi\\ in one rule, and
every category but the reference enters the same way.

## Notation

\\K\\ is the number of categories, \\\eta \in \mathbb{R}^{K-1}\\ the
free vector, \\\pi\\ the probability vector and \\\delta\_{ab}\\ the
Kronecker delta. The reference category is the last, \\K\\.

## References

Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*.
Chapman and Hall, London. The additive log-ratio chart is chapter 6.

## See also

[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
which is one of these per row,
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
for the map and its inverse, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md),
whose seven-check battery is what a non-matrix family gets.

## Examples

``` r
# Four categories, three free values.
s <- simplex(4)
s@free_names
#> [1] "alr1" "alr2" "alr3"
eta <- c(0.5, -0.2, 1.1)
pi <- param_value(s, eta)
round(pi, 5)
#>      p1      p2      p3      p4 
#> 0.25476 0.12651 0.46421 0.15452 
sum(pi)
#> [1] 1

# It is the closed form, with the last category as reference.
all.equal(unname(pi),
          c(exp(eta), 1) / (1 + sum(exp(eta))))
#> [1] TRUE

# The round trip closes exactly.
max(abs(param_free(s, pi) - eta))
#> [1] 2.220446e-16

# The first derivative is the covariance of a categorical indicator.
d1 <- param_d1(s, eta)
round(d1[["alr1"]], 6)
#> [1]  0.189858 -0.032230 -0.118262 -0.039366
round(pi * ((seq_along(pi) == 1) - pi[1]), 6)
#>        p1        p2        p3        p4 
#>  0.189858 -0.032230 -0.118262 -0.039366 

# Every derivative component sums to zero over the value index, at every
# order, because sum(pi) is the constant 1.
vapply(1:4, function(k) {
  dk <- do.call(paste0("param_d", k), list(s, eta))
  max(abs(vapply(dk, sum, numeric(1))))
}, numeric(1))
#> [1] 6.938894e-18 6.938894e-18 1.040834e-17 1.734723e-17

# A large free value saturates instead of overflowing, and the vector still
# sums to exactly 1.
v <- param_value(s, c(800, 0, 0))
c(p1 = v[[1]], pK = v[[4]], total = sum(v))
#>    p1    pK total 
#>     1     0     1 
```
