# Second Derivatives of a Parameter's Value

Differentiates the map twice, returning the \\d(d+1)/2\\ distinct
components \\\partial^2 V / \partial \eta_k \partial \eta_l\\, each
shaped like the value. A mixed partial does not depend on the order of
differentiation, so only the unordered pairs are returned; the caller
reads the \\(l, k)\\ entry off the \\(k, l)\\ one. A Newton step in the
free vector needs these, and so does the observed information of a model
whose covariance is parametrized this way.

## Usage

``` r
param_d2(s, eta, ...)
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

A named list of `choose(s@n_free + 1, 2)` entries, keyed as
`param_tuple_names(s)` and in that order. Each entry has the shape of
the value, symmetric for a matrix family.

## Keys and their order

The list is keyed by
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md),
which puts the \\d\\ diagonal pairs first and then the \\d(d-1)/2\\
off-diagonal ones in lexicographic order. Diagonal first is what a
consumer filling a Hessian wants, and the ordering is part of the
interface:
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
returns the index pairs in exactly the same order, so a caller can walk
the two together.

Both come from one enumeration, and neither is produced by taking a key
apart. Splitting `"log_L1:log_L2"` on `":"` looks equivalent and is not:
a free value whose own label contains the separator splits into the
wrong number of pieces, and the failure is silent.

## Exact, or one stencil

Every family here is closed form. The fallback on
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
takes one difference and no more, and which quantity it differences
depends on what the family already supplies. Where
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
is analytic it differences that once, in the other index. Where
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
is itself numerical it works on
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
directly: a three-point second difference where the two indices
coincide, and one difference in each of the two components where they
differ. A mixed derivative in two different variables is one stencil
however it is written; the nesting the toolkit forbids is two
differences in the same variable.

## Notation

\\\eta\\ is the free vector, of length \\d = \\ `s@n_free`, and \\V\\
the value
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns.

## See also

[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
for the first order,
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for the higher ones,
[`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
and
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
for the keys, and
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
for the second derivative of the log-determinant, which is not the trace
of this.

## Examples

``` r
s <- log_cholesky(2)
eta <- c(0.2, -0.1, 0.4)
d2 <- param_d2(s, eta)
names(d2)
#> [1] "log_L1:log_L1" "log_L2:log_L2" "L2.1:L2.1"     "log_L1:log_L2"
#> [5] "log_L1:L2.1"   "log_L2:L2.1"  

# The diagonal pairs come first, and the keys match the index tuples.
param_tuple_indices(s, 2)[1:3]
#> [[1]]
#> [1] 1 1
#> 
#> [[2]]
#> [1] 2 2
#> 
#> [[3]]
#> [1] 3 3
#> 

# A second difference of the map agrees with the L2.1 diagonal component.
h <- 1e-4
e <- c(0, 0, h)
fd <- (param_value(s, eta + e) - 2 * param_value(s, eta) +
         param_value(s, eta - e)) / h^2
max(abs(d2[["L2.1:L2.1"]] - fd))
#> [1] 1.004952e-08

# log_cholesky is quadratic in a below-diagonal free value, so the second
# derivative there is constant and the third vanishes.
d2[["L2.1:L2.1"]]
#>    v1 v2
#> v1  0  0
#> v2  0  2
max(abs(param_d3(s, eta)[["L2.1:L2.1:L2.1"]]))
#> [1] 0
```
