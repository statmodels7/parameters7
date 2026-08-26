# Fourth Derivatives of a Parameter's Value

Differentiates the map four times, returning the distinct components

\$\$\frac{\partial^{4} V(\eta)}
{\partial\eta_k\\\partial\eta_l\\\partial\eta_m\\\partial\eta_n},\$\$

each shaped like the value. Fourth order is where the contract stops. A
fourth-order chain rule through a link needs exactly this much, and
nothing in the toolkit asks for a fifth.

## Usage

``` r
param_d4(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md).

- eta:

  A numeric vector of length `s@n_free`, finite in every entry.

- ...:

  Passed to the method. No method in this package reads it.

## Value

A named list of `choose(s@n_free + 3, 4)` entries, keyed as
`param_tuple_names(s, 4)` and in that order. Each entry has the shape of
the value.

## Keys

The derivative is symmetric in its indices, so the \\\binom{d+3}{4}\\
unordered quadruples are returned, keyed by `param_tuple_names(s, 4)` in
lexicographic order, with `param_tuple_indices(s, 4)` giving the index
quadruples in the same order. At \\d = 3\\ there are fifteen, at \\d =
6\\ a hundred and twenty-six, so the list grows quickly in the size of
the matrix.

## Exact, or one stencil

Every family here is closed form. The fallback on
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
applies one product stencil per quadruple, with a factor per distinct
index of the width its multiplicity calls for, evaluated in one pass
over
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md).
At fourth order the rounding of a difference grows as \\\varepsilon /
h^4\\, so a numerical answer here is the least accurate of the four
orders; a family fitted in earnest is better served by writing the
closed form out.

## Notation

\\\eta\\ is the free vector, of length \\d = \\ `s@n_free`, and \\V\\
the value
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns.

## See also

[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
for the order below,
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
for the keys, and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for the fourth derivative of the log-determinant.

## Examples

``` r
# For a scalar matrix, exp(eta) times the identity, every order is the
# matrix again.
s <- scalar_matrix(2)
all.equal(param_d4(s, 0.3)[[1]], param_value(s, 0.3))
#> [1] TRUE

# The list is over unordered quadruples, so it is shorter than d^4.
q <- log_cholesky(3)
c(returned = length(param_d4(q, rep(0.1, 6))),
  unordered = choose(6 + 3, 4), ordered = 6^4)
#>  returned unordered   ordered 
#>       126       126      1296 

# log_cholesky is quadratic in a below-diagonal free value, so any quadruple
# repeating one of those three or more times is exactly zero.
d4 <- param_d4(log_cholesky(2), c(0.2, -0.1, 0.4))
max(abs(d4[["L2.1:L2.1:L2.1:L2.1"]]))
#> [1] 0
max(abs(d4[["log_L1:L2.1:L2.1:L2.1"]]))
#> [1] 0
```
