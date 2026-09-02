# Construct a Diagonal Parameter

Returns an object holding the map \\M = \mathrm{diag}(h(\eta_1), \dots,
h(\eta_p))\\, a diagonal matrix whose \\p\\ entries are positive, each
carried from one free value by a scalar link. Use it for a covariance
that assumes independence, which is \\p\\ free values against
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)'s
\\p(p+1)/2\\, and for a set of independent variance components.

The choice of link is yours, and it is a real choice: the derivatives
the family reports are the link's own, so the curvature an optimizer
sees on the free scale changes with it.
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
is the same object with one free value shared by every entry.

## Usage

``` r
diagonal_matrix(dimension, link = linkfunctions7::log_link())
```

## Arguments

- dimension:

  The side \\p\\ of the matrix. A single positive whole number, finite
  and at least 1; anything else throws
  `'dimension' must be a single positive integer.`

- link:

  A linkfunctions7 link carrying the free scale onto the positive
  entries,
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default. It must map **onto** the positive half line, so
  `identity_link()` is rejected with its bounds in the message, and
  **from** the whole real line, since the free vector is unconstrained
  by design.
  [`check_positive_link()`](https://statmodels7.github.io/parameters7/reference/check_positive_link.md)
  enforces both at construction, reading the link's own lower bound for
  the first and
  [`linkfunctions7::eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.html)
  for the second. `softplus_link()` satisfies both, and so does a
  bounded link such as `logit_link()`, whose range \\(0, 1)\\ is
  positive and which gives a positive definite matrix with entries
  below 1. `sqrt_link()`, `inverse_link()`, `inverse_sq_link()` and
  `power_link()` at a positive exponent satisfy only the first and are
  rejected: their predictor scale is \\(0, \infty)\\, so the map is even
  in \\\eta\\ and the round trip returns \\\lvert \eta \rvert\\.

## Value

An object of class
[`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md),
with `n_free` equal to `dimension`, `free_names` tagged by the link,
`rank` equal to `dimension`, an empty `null_basis`, `param_name`
`"diag"`, and `param_params` holding `link` and `shared = FALSE`.

## Where the two packages meet

The Jacobian of a diagonal map is diagonal, which is exactly the
contract a scalar link satisfies, so the linkfunctions7 objects are
reused as they are and their exact derivatives to fourth order come with
them. Nothing here rederives a chain rule.

## The free names record the link

`free_names` is `log_d1`, `log_d2`, ... under the default, and
`softplus_d1`, `softplus_d2`, ... under a softplus link. A label names
the coordinate, never the quantity it produces, so a consumer that
flattens the free vector into scalars with identity links reports a
number on the scale it is really on.

## The log-determinant, and when it is linear

\$\$\log\|M\| = \sum\_{i=1}^{p} \log h(\eta_i),\$\$

a sum of functions of one free value each, so every mixed derivative of
it is exactly zero at every order. Under the **log** link \\\log h(\eta)
= \eta\\, so the gradient is a vector of ones and the second, third and
fourth derivatives all vanish. Under any other link they do not: a
square-root link at \\\eta = (1, 2)\\ gives a second derivative of
\\(-2, -0.5)\\ on the diagonal.

## Notation

\\\eta\\ is the free vector, of length \\d = p\\, and \\p\\ the side of
the matrix. \\h = g^{-1}\\ is the inverse link, carrying a free value
onto a positive diagonal entry, and \\h'\\, \\h''\\ its derivatives in
\\\eta\\.

## See also

[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
for one shared value,
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
for an unstructured matrix,
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
to give a correlation its own diagonal scale, and
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
for the map.

## Examples

``` r
# Three independent variances, on the log scale by default.
s <- diagonal_matrix(3)
s@free_names
#> [1] "log_d1" "log_d2" "log_d3"
round(param_value(s, c(0, 0.5, -0.5)), 4)
#>    v1     v2     v3
#> v1  1 0.0000 0.0000
#> v2  0 1.6487 0.0000
#> v3  0 0.0000 0.6065

# The entries are exp() of the free values, and the round trip closes.
all.equal(unname(diag(param_value(s, c(0, 0.5, -0.5)))),
          exp(c(0, 0.5, -0.5)))
#> [1] TRUE
max(abs(param_free(s, param_value(s, c(0.2, -0.3, 0.7))) -
        c(0.2, -0.3, 0.7)))
#> [1] 1.110223e-16

# Any link carrying the whole free line onto positive entries serves, and
# the free names say which one was used.
r <- diagonal_matrix(2, link = linkfunctions7::softplus_link())
r@free_names
#> [1] "softplus_d1" "softplus_d2"
round(param_value(r, c(1, 2)), 4)
#>        v1     v2
#> v1 1.3133 0.0000
#> v2 0.0000 2.1269

# A link onto the whole line is refused, an entry having to be positive,
# and so is one defined on part of the line, the map being even there.
try(diagonal_matrix(2, link = linkfunctions7::identity_link()))
#> Error : 'link' maps onto (-Inf, Inf), which is not inside the positive half
#>   line. A diagonal entry of a positive definite matrix is positive.
try(diagonal_matrix(2, link = linkfunctions7::sqrt_link()))
#> Error : 'link' is defined on the predictors (0, Inf) and not on the whole
#>   real line. The free vector is unconstrained by design, so a map
#>   reaching only part of the line gives the same matrix at two free
#>   vectors and the round trip returns |eta|.
#>   log_link() and softplus_link() are the shipped links onto the
#>   positive half line from the whole of it.

# The log-determinant is a sum over the entries, and under the log link it
# is linear, so its higher derivatives vanish exactly.
eta <- c(0, 0.5, -0.5)
all.equal(param_logdet(s, eta), sum(eta))
#> [1] TRUE
param_dlogdet(s, eta)
#> log_d1 log_d2 log_d3 
#>      1      1      1 
max(abs(param_d2logdet(s, eta)))
#> [1] 0

# Under another link they do not.
param_d2logdet(r, c(1, 2))
#> softplus_d1:softplus_d1 softplus_d2:softplus_d2 softplus_d1:softplus_d2 
#>              -0.1601732              -0.1221289               0.0000000 
```
