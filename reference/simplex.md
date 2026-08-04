# Construct a Simplex Parameter

A probability vector of length \\K\\ on the open simplex, carried by the
additive log-ratio map with the last category as reference: \$\$\pi_a =
\dfrac{e^{\eta_a}}{1 + \sum_b e^{\eta_b}}, \qquad \pi_K = \dfrac{1}{1 +
\sum_b e^{\eta_b}},\$\$ with \\\eta \in \mathbb{R}^{K-1}\\
unconstrained.

## Usage

``` r
simplex(n_cat)
```

## Arguments

- n_cat:

  The number of categories \\K\\, at least 2.

## Value

An object of class
[`SimplexParam`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md).

## Details

The derivatives of the map close over the value itself. With \\Z = 1 +
\sum_b e^{\eta_b}\\, each \\\pi_a\\ with \\a \< K\\ is a derivative of
\\\log Z\\, so the derivative tensors of \\\pi\\ are the cumulants of a
categorical indicator and follow from the one rule \\\partial_b \pi_a =
\pi_a(\delta\_{ab} - \pi_b)\\ applied repeatedly. All four orders are
closed form.

The inverse is \\\eta_a = \log(\pi_a / \pi_K)\\, exact; a vector outside
the open simplex, or one that does not sum to one, is refused rather
than repaired, because a silent renormalisation would mask the caller's
defect.

Stick-breaking is the other common chart and was not chosen: its
derivatives chain through \\K - 1\\ nested logistic maps and are not
symmetric in the categories, while the additive log-ratio's close over
\\\pi\\ in one rule.

## See also

[`transition_matrix`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
[`param_value`](https://statmodels7.github.io/parameters7/reference/param_value.md)

## Examples

``` r
s <- simplex(3)
round(param_value(s, c(0.5, -0.2)), 4)
#>     p1     p2     p3 
#> 0.4755 0.2361 0.2884 

# the round trip closes exactly
eta <- c(0.5, -0.2)
max(abs(param_free(s, param_value(s, eta)) - eta))
#> [1] 2.775558e-17
```
