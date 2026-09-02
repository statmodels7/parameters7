# Construct a Transition Matrix Parameter

Returns an object holding a \\K \times K\\ row-stochastic matrix, the
transition matrix of a Markov chain on \\K\\ states, with each row an
independent
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
in the additive log-ratio chart. \\K(K-1)\\ free values in all, and
every free vector gives positive entries with rows summing to exactly 1,
so a chain can be estimated without a constraint.

## Usage

``` r
transition_matrix(n_state)
```

## Arguments

- n_state:

  The number of states \\K\\, **at least 2**. A single integer; `1`, a
  fraction, `NA` and a vector all throw
  `'n_state' must be a single integer of at least 2.` A one-state chain
  has nothing to estimate.

## Value

An object of class
[`TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md),
with `n_free` equal to \\K(K-1)\\, `free_names` `alr1.1`, `alr1.2`, ...,
`alr{K}.{K-1}` row by row, `param_name` `"transition_matrix"`, and
`param_params` holding `n_state`. No `dimension`, `rank` or
`null_basis`.

## Rows, not columns

The rows live on the simplex, so row \\i\\ holds the probabilities of
moving **from** state \\i\\. A transition matrix in this convention acts
on row vectors of probabilities, \\\pi\_{t+1} = \pi_t P\\; the column
convention is its transpose.

## The rows are independent

Each row is parametrized by its own \\K-1\\ free values and no others,
so every derivative array is **block diagonal by row**: a component
pairing free values of two different rows is exactly the zero matrix,
and a component inside row \\i\\ is zero everywhere outside row \\i\\.
The implementation evaluates the simplex kernels row by row and never
stores those zeros.

Measured at \\K = 3\\: the second-derivative component `alr1.1:alr2.1`,
spanning two rows, is 0 exactly, while `alr1.1:alr1.1` is not.

## Reading the free vector

It runs row by row, and the names `alr{i}.{j}` say which row and which
chart coordinate, the row.column convention
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
uses. The `alr` records the chart, so a free value of 0.5 is not a
probability of 0.5: it is a log ratio against the row's last state.

## What it has and has not got

All four derivative orders are closed form, inherited from
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)'s
cumulant recursion. There is no log-determinant, no solve and no factor:
the value is not symmetric, and the object is a
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
and never a
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
so those generics have no method and fail at dispatch.
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
runs its seven-check battery instead of the nine-check one, and the
identity it tests here is that every derivative component sums to zero
along each row.

## Notation

\\K\\ is the number of states, \\\eta \in \mathbb{R}^{K(K-1)}\\ the free
vector and \\P\\ the transition matrix, with \\P\_{ij}\\ the probability
of moving from state \\i\\ to state \\j\\. The reference category of
each row is its last state.

## See also

[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md),
which is one row of this and carries the derivative recursion,
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
for the map and its inverse, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
for the battery a non-matrix family gets.

## Examples

``` r
# Three states: six free values, two per row.
s <- transition_matrix(3)
c(n_free = s@n_free)
#> n_free 
#>      6 
s@free_names
#> [1] "alr1.1" "alr1.2" "alr2.1" "alr2.2" "alr3.1" "alr3.2"

set.seed(4)
eta <- rnorm(s@n_free)
m <- param_value(s, eta)
round(m, 4)
#>        s1     s2     s3
#> s1 0.4399 0.2059 0.3542
#> s2 0.4641 0.3455 0.1904
#> s3 0.6317 0.2452 0.1231

# Every row is a probability distribution, exactly.
rowSums(m)
#> s1 s2 s3 
#>  1  1  1 

# The round trip closes exactly.
max(abs(param_free(s, m) - eta))
#> [1] 2.220446e-16

# The rows are independent, so a derivative in a row-1 free value is zero
# everywhere outside row 1.
round(param_d1(s, eta)[["alr1.1"]], 4)
#>        s1      s2      s3
#> s1 0.2464 -0.0906 -0.1558
#> s2 0.0000  0.0000  0.0000
#> s3 0.0000  0.0000  0.0000

# And a second derivative spanning two rows is exactly the zero matrix.
d2 <- param_d2(s, eta)
c(across_rows = max(abs(d2[["alr1.1:alr2.1"]])),
  within_row = max(abs(d2[["alr1.1:alr1.1"]])))
#> across_rows  within_row 
#>  0.00000000  0.02960676 

# There is no log-determinant to ask for: the value is not symmetric.
try(param_logdet(s, eta))
#> Error : Can't find method for `param_logdet(<parameters7::TransitionMatrixParam>)`.
```
