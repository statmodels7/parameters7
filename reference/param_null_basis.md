# An Orthonormal Basis of a Null Space, and the Rank That Goes With It

Returns the rank and an orthonormal basis of the common null space of
one or more symmetric positive semidefinite matrices. Each matrix is
scaled to unit maximum entry, the scaled matrices are stacked into one
tall matrix, and the rank and the null space are read off its singular
value decomposition. Every
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
calls it once at construction to fill its `rank` and `null_basis`
properties.

## Usage

``` r
param_null_basis(mats, tol = 1e-10)
```

## Arguments

- mats:

  A list of symmetric numeric matrices, all of the same side, or a
  single matrix, which is wrapped in a list. An empty list throws
  `'mats' must not be empty.`; the matrices themselves are not checked
  for symmetry or definiteness, the callers being the package's own
  constructors.

- tol:

  The relative tolerance below which a singular value counts as zero: a
  singular value \\s_j\\ is null when \\s_j \le \mathrm{tol} \cdot
  \max_i s_i\\. Defaults to `1e-10`. The default sits well below the
  singular values a genuine rank carries and well above the \\10^{-16}\\
  the zero directions of a stacked, normalized set of matrices reach, so
  the gap between the two is wide and the exact value is not delicate.

## Value

A list with two components

- `rank`:

  an integer, the number of singular values above the tolerance.

- `null_basis`:

  a `p` by `p - rank` numeric matrix with orthonormal columns, spanning
  the common null space. It has zero columns when the stack has full
  rank.

## Reading the rank off the components

The null space of a sum of positive semidefinite matrices is the
intersection of their null spaces, and the intersection is the null
space of the stack, so reading the answer off the stack is exact.
Reading it off an assembled combination \\\lambda_1 P_1 + \lambda_2
P_2\\ is not, because a count of small eigenvalues is not scale
invariant: a component whose weight is small contributes eigenvalues
below the tolerance, and they are then counted as null directions that
are not there.

Measured on the tensor-product penalty of two second-difference
penalties over 4 and 8 coefficients, whose true rank is 28 out of 32.
Counting the eigenvalues of the assembled sum above a relative tolerance
of \\10^{-10}\\ answers 28 while the two weights are within \\10^{8}\\
of each other, and 24 once the ratio reaches \\10^{10}\\: four
directions the penalty does penalize are read as null. Smoothing
parameters ten orders of magnitude apart are an ordinary fitted model,
not a pathology. The stacked route answers 28 at every ratio, and the
basis it returns is annihilated by the assembled matrix to \\3 \times
10^{-16}\\ relative even at the worst one.

## The scaling

Each component is divided by its largest absolute entry before stacking,
so that a component contributed in different units still gets an equal
vote. A component that is identically zero is passed through as zero and
contributes nothing, which is correct: its null space is the whole
space.

## Cost

One singular value decomposition of an \\(kp) \times p\\ matrix for
\\k\\ components of side \\p\\, so \\O(k p^3)\\. It runs once per
constructed parameter and never inside a fit.

## Notation

\\p\\ is the side of the matrices and \\k\\ their number. The **null
basis** is a matrix whose columns span \\\\v : M v = 0 \text{ for every
} M\\\\, and the **rank** is \\p\\ minus the dimension of that space.

## See also

[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
whose `rank` and `null_basis` properties this fills, and the two
families that can be rank deficient and so are the ones whose answer is
not trivial:
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md).

## Examples

``` r
# A second-difference penalty on six coefficients has rank 4, its null space
# being the constants and the straight lines.
d <- diff(diag(6), differences = 2)
nb <- param_null_basis(crossprod(d))
nb$rank
#> [1] 4
dim(nb$null_basis)
#> [1] 6 2

# The claim, checked: a constant and a line are in the span and a quadratic
# is not. Projecting onto the basis leaves the first two untouched.
proj <- function(v) nb$null_basis %*% (t(nb$null_basis) %*% v)
vapply(list(rep(1, 6), 1:6, (1:6)^2),
       function(v) sqrt(sum((v - proj(v))^2)) / sqrt(sum(v^2)), numeric(1))
#> [1] 1.351405e-15 6.105395e-16 1.281025e-01

# Why the components are taken separately. Two marginal penalties on a
# 4 x 8 tensor product: the pair has rank 28 out of 32.
P <- function(m) crossprod(diff(diag(m), differences = 2))
P1 <- kronecker(P(4), diag(8))
P2 <- kronecker(diag(4), P(8))
param_null_basis(list(P1, P2))$rank
#> [1] 28

# Counting small eigenvalues of the assembled sum agrees while the two
# weights are comparable and loses four directions when they are not.
count <- function(M, tol = 1e-10) {
  e <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  sum(e > tol * max(e))
}
c(equal = count(P1 + P2), ratio_1e10 = count(P1 + 1e10 * P2))
#>      equal ratio_1e10 
#>         28         24 

# The stacked answer stays correct there, and its basis is annihilated.
nb2 <- param_null_basis(list(P1, P2))
M <- P1 + 1e10 * P2
max(abs(M %*% nb2$null_basis)) / max(abs(M))
#> [1] 6.717643e-16
```
