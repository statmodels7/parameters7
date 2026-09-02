# Construct a Block Diagonal of Matrix Parameters

Returns an object holding \\M(\eta) = \mathrm{diag}(S_1(\eta_1), \ldots,
S_B(\eta_B))\\ for matrix parameters \\S_1, \ldots, S_B\\: a
block-diagonal matrix whose blocks are **distinct families**, each
governed by its own stretch of the free vector.

It is the covariance (or precision) of several independent groups of
coefficients whose structures differ, as a model carrying more than one
random-effect term has.

## Usage

``` r
block_diag(...)
```

## Arguments

- ...:

  One or more objects inheriting from
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  or a single list of them. Named arguments supply the block labels,
  which must be unique. A block that is not a `matrix_parameter` is
  rejected by position, so `block_diag(ar1(3), simplex(3))` reports that
  block 2 does not inherit from it.

## Value

An object of class
[`BlockDiagParam()`](https://statmodels7.github.io/parameters7/reference/BlockDiagParam.md),
with `dimension`, `n_free` and `rank` the sums of the blocks',
`free_names` the blocks' own prefixed by the labels, and `null_basis`
the block diagonal of the blocks'.

## Nothing is rederived, and the reason

The free values of one block do not enter another, so

\$\$\partial_k M = \mathrm{diag}(0, \ldots, \partial_k S_b, \ldots, 0),
\qquad \log\lvert M \rvert\_{+} = \sum_b \log\lvert S_b \rvert\_{+},
\qquad M^{-1} = \mathrm{diag}(S_1^{-1}, \ldots, S_B^{-1}),\$\$

and a derivative whose indices do not all belong to one block is
**identically zero**, at every order and for the log-determinant as well
as for the value. Measured on `block_diag(log_cholesky(2), ar1(3))`,
whose five free values split 3 and 2: 6 of the 15 second-order
components are cross-block, 21 of 35 at third order and 50 of 70 at
fourth, and every one of them is exactly 0. What is left is fetched from
the block and placed in the rows and columns that block occupies.

The rank is the sum of the blocks' ranks and the null basis is their
block diagonal, both read from the components and never from an
assembled matrix, which is the rule this package follows everywhere: a
rank is a property of the family, and counting eigenvalues of the
assembled matrix would make it a property of the arithmetic. A deficient
block is therefore admitted, and the composite reports the deficiency:
with a rank-one \\P\\ of side 2, `block_diag(scaled_matrix(P), ar1(3))`
has rank 4 of 5 and a null basis of one column,
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
returns the log pseudo-determinant, and
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
is refused by the generic.

## Against kron_identity()

[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
repeats **one** block \\m\\ times and they share a single free vector,
so its `n_free` does not grow with \\m\\. Here the blocks are different
objects and their free vectors are concatenated, so the composite has
\\\sum_b d_b\\ free values. Use this one where the groups have different
structures, and that one where they have the same structure and the same
parameters.

## Labels

The free names are prefixed by the block's label, since two blocks of
the same family would otherwise report the same names and the class
requires them to be unique. Labels come from the names of the arguments
where they are given, and are `b1`, `b2`, ... otherwise. `param_name`
records the blocks' own names, as in `"blockdiag(log_cholesky, ar1)"`.

## Notation

\\B\\ is the number of blocks, \\S_b\\ the \\b\\-th block's map, \\d_b\\
its number of free values, and \\\eta_b\\ the stretch of the free vector
it owns.

## See also

[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
for identical blocks,
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
for the other two compositions, and
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the contract every block meets.

## Examples

``` r
s <- block_diag(subject = log_cholesky(2), time = ar1(3))
c(dimension = s@dimension, n_free = s@n_free)
#> dimension    n_free 
#>         5         5 
s@free_names
#> [1] "subject_log_L1" "subject_log_L2" "subject_L2.1"   "time_log_scale"
#> [5] "time_z_rho"    

eta <- c(0.1, -0.2, 0.3, log(2), atanh(0.5))
M <- param_value(s, eta)

# Each block is its own family's value, and the off-diagonal is exactly zero.
c(subject = max(abs(M[1:2, 1:2] - param_value(log_cholesky(2), eta[1:3]))),
  time = max(abs(M[3:5, 3:5] - param_value(ar1(3), eta[4:5]))),
  off = max(abs(M[1:2, 3:5])))
#> subject    time     off 
#>       0       0       0 

# The log-determinant is the sum of the blocks'.
c(composite = param_logdet(s, eta),
  sum_of_blocks = param_logdet(log_cholesky(2), eta[1:3]) +
                  param_logdet(ar1(3), eta[4:5]))
#>     composite sum_of_blocks 
#>      1.304077      1.304077 

# And 50 of the 70 fourth-order components are exactly zero, their indices
# not all belonging to one block.
owner <- s@param_params$owner
cross <- vapply(param_tuple_indices(s, 4L),
                function(t) length(unique(owner[t])) > 1L, logical(1))
c(cross_block = sum(cross), of = length(cross),
  largest = max(abs(unlist(param_d4(s, eta)[cross]))))
#> cross_block          of     largest 
#>          50          70           0 

# The round trip closes.
max(abs(param_free(s, M) - eta))
#> [1] 6.938894e-17
```
