# Names of the Distinct Derivative Components

Returns the keys of a derivative list: one per unordered tuple of free
values, built by pasting the free names together with `":"`. At order 2
these are the keys of
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
and
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
with the diagonal pairs first; at orders 3 and 4 the keys of
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md),
in lexicographic order. Use it to look a component up by name, or to
walk a derivative list beside the index tuples
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
returns.

## Usage

``` r
param_tuple_names(s, order = 2L)
```

## Arguments

- s:

  An object inheriting from class
  [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md),
  whose `free_names` supply the labels and whose `n_free` supplies
  \\d\\.

- order:

  The derivative order. `2` by default, which is the order
  [`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
  and
  [`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
  use; `3` and `4` name the components of
  [`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
  and
  [`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md).
  `1` is accepted and returns the free names themselves. Anything else
  throws `'order' must be 1, 2, 3 or 4.`, the contract stopping at
  fourth order.

## Value

A character vector of `choose(s@n_free + order - 1, order)` names, in
the same order as the list they key. At order 2 the `s@n_free` diagonal
pairs come first.

## How many there are

A mixed partial does not depend on the order of differentiation, so the
components at order \\k\\ are the multisets of size \\k\\ drawn from the
\\d\\ free values, of which there are

\$\$\binom{d + k - 1}{k},\$\$

against \\d^k\\ ordered tuples. That is the length of the vector
returned and of every derivative list of that order. For \\d = 6\\,
which is a \\3 \times 3\\ unstructured covariance, order 4 has 126
components against 1296 ordered ones.

## Why the keys and the tuples come from one enumeration

Nothing in the toolkit recovers an index by splitting a key apart.
Taking `"log_L1:log_L2"` and splitting on `":"` works until a free name
contains the separator itself, and then it yields the wrong number of
pieces and the failure is silent. Generating the names and the indices
from one enumeration cannot be fooled that way. `param_tuple_names()`
calls
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
and labels what it gets, so the two agree by construction.

## Notation

\\d\\ is the length of the free vector, `s@n_free`, and \\k\\ the
derivative order.

## See also

[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
for the index tuples in the same order, and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md),
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md),
[`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
[`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
and
[`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
for the lists these key.

## Examples

``` r
s <- log_cholesky(2)
param_tuple_names(s)
#> [1] "log_L1:log_L1" "log_L2:log_L2" "L2.1:L2.1"     "log_L1:log_L2"
#> [5] "log_L1:L2.1"   "log_L2:L2.1"  
param_tuple_names(s, 3)
#>  [1] "log_L1:log_L1:log_L1" "log_L1:log_L1:log_L2" "log_L1:log_L1:L2.1"  
#>  [4] "log_L1:log_L2:log_L2" "log_L1:log_L2:L2.1"   "log_L1:L2.1:L2.1"    
#>  [7] "log_L2:log_L2:log_L2" "log_L2:log_L2:L2.1"   "log_L2:L2.1:L2.1"    
#> [10] "L2.1:L2.1:L2.1"      

# They really are the keys of the derivative list, in order.
identical(param_tuple_names(s, 4), names(param_d4(s, c(0.2, -0.1, 0.4))))
#> [1] TRUE

# The count is over unordered tuples.
q <- log_cholesky(3)
vapply(1:4, function(k) length(param_tuple_names(q, k)), integer(1))
#> [1]   6  21  56 126
choose(6 + 1:4 - 1, 1:4)
#> [1]   6  21  56 126

# Order 2 puts the diagonal pairs first, which is how a consumer filling a
# Hessian wants them.
param_tuple_names(s, 2)
#> [1] "log_L1:log_L1" "log_L2:log_L2" "L2.1:L2.1"     "log_L1:log_L2"
#> [5] "log_L1:L2.1"   "log_L2:L2.1"  
```
