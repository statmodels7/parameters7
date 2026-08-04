# Construct a Transition Matrix Parameter

A \\K \times K\\ row-stochastic matrix – the transition matrix of a
Markov chain on \\K\\ states – with each row an independent
[`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.md)
in the additive log-ratio chart, so \\K(K-1)\\ free values in all.

## Usage

``` r
transition_matrix(n_state)
```

## Arguments

- n_state:

  The number of states \\K\\, at least 2.

## Value

An object of class
[`TransitionMatrixParam`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md).

## Details

The rows are independent in the parametrisation, so every derivative
tensor is block diagonal by row: a component pairing free values of two
different rows is exactly zero, and the implementation evaluates the
simplex kernels row by row rather than storing those zeros.

The free vector runs row by row, and the names `alr{i}.{j}` say which
row and which chart coordinate, the row.column convention the
log-Cholesky names already use. Rows, not columns, live on the simplex:
a transition matrix acts on row vectors of probabilities, and the column
convention is its transpose.

## See also

[`simplex`](https://statmodels7.github.io/parameters7/reference/simplex.md)

## Examples

``` r
s <- transition_matrix(3)
eta <- rnorm(s@n_free)
rowSums(param_value(s, eta))
#> s1 s2 s3 
#>  1  1  1 
```
