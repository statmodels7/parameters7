# The Abstract Class of a Symmetric Matrix Parameter

Extends
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
to the branch whose value is a symmetric positive semidefinite \\p
\times p\\ matrix. Beyond the map itself, such a family can answer the
four things a Gaussian likelihood asks of a covariance or a precision:
the log-determinant, a solve against a right-hand side, a factor, and
the rank and null space when the matrix is singular. Those four generics
are registered here, so a subclass supplying only
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
inherits all of them.

The class is abstract. Every matrix family in the package returns a
subclass of it:
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
and the four compositions.

## Usage

``` r
matrix_parameter(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list(),
  dimension = integer(0),
  rank = integer(0),
  null_basis = integer(0),
  role = character(0)
)
```

## Arguments

- param_name:

  A single character string naming the family.

- n_free:

  The length \\d\\ of the free vector: a single non-negative integer,
  agreeing with `length(free_names)`.

- free_names:

  A character vector of length `n_free`, one label per free value, in
  the order the free vector holds them. Must be unique.

- param_params:

  A list of whatever the family needs in order to evaluate itself, read
  only by that family's own methods.

- dimension:

  The side \\p\\ of the matrix: a single integer, no smaller than 1. The
  validator rejects a vector and a value below 1.

- rank:

  The rank of the matrix the family produces, a single integer in
  `0:dimension`. It is a property of the family, so a family whose value
  is positive definite at every \\\eta\\ declares \\p\\ here.

- null_basis:

  A `dimension` by `dimension - rank` numeric matrix whose columns are
  an orthonormal basis of the common null space. Use
  [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
  to obtain one, or `matrix(numeric(0), dimension, 0)` for a full-rank
  family. The validator rejects any other shape, and reports both the
  rank and the shape when the two disagree.

- role:

  A single string, one of `"covariance"`, `"precision"` or `"either"`,
  recording which side of a model the matrix parametrizes. **No numeric
  result depends on it.** It is carried because the family name does not
  record it: the same
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  serves either side, and a consumer that prefixes a free name with the
  matrix it describes needs to know which.
  [`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
  reads it to give a composite the common role of its blocks, or
  `"either"` when they disagree, and
  [`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
  copies it.

## Value

An object of class `matrix_parameter`, which is a
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
with four further properties

- `dimension`:

  integer, the side \\p\\.

- `rank`:

  integer in `0:dimension`.

- `null_basis`:

  a `dimension` by `dimension - rank` matrix with orthonormal columns.

- `role`:

  character, as supplied.

plus the four it inherits, `param_name`, `n_free`, `free_names` and
`param_params`. The class is abstract, so a useful object comes from a
constructor such as
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md).

## Rank and null space are fixed at construction

They are properties of the family, the same at every point \\\eta\\.
Scaling a matrix by a positive number leaves its null space alone, and
the null space of a sum of positive semidefinite matrices is the
intersection of theirs, so no free value moves it. Recording them once
at construction is exact, and it is also the only stable way to obtain
them: counting the small eigenvalues of an assembled matrix is not scale
invariant, and reads a component whose weight is small as a null
direction.
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
computes them from the components and carries the measurement.

Most families here are full rank, so `rank` is \\p\\ and `null_basis`
has zero columns.
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
admits a deficient fixed matrix, and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
a deficient sum, and both then report the deficiency.

## What a non-matrix family does instead

[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
and
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
inherit
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
directly. They hold no `dimension`, `rank` or `null_basis`, and
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
has no method for them, so asking for the log-determinant of a
probability vector fails at dispatch. The absence is structural.

## Notation

\\\eta\\ is the free vector, the point on the unconstrained scale, and
\\d\\ its length. \\p\\ is the side of the matrix. \\M\\ is the matrix
the map produces, \\\Sigma\\ when it is read as a covariance and
\\\Omega\\ when it is read as a precision.

## See also

[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
for the base class and the naming convention for free values.
[`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md),
[`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
and
[`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
for the four quantities this branch adds.

## Examples

``` r
# A matrix family is both a `parameter` and a `matrix_parameter`.
s <- log_cholesky(3)
c(parameter = S7::S7_inherits(s, parameter),
  matrix_parameter = S7::S7_inherits(s, matrix_parameter))
#>        parameter matrix_parameter 
#>             TRUE             TRUE 

# The four properties this branch adds.
c(dimension = s@dimension, rank = s@rank, null_columns = ncol(s@null_basis))
#>    dimension         rank null_columns 
#>            3            3            0 
s@role
#> [1] "either"

# Belonging to this branch is what gives the family a log-determinant, and
# it agrees with the eigenvalues of the matrix itself.
eta <- c(0.1, -0.2, 0.3, 0.5, -0.4, 0.2)
all.equal(param_logdet(s, eta),
          sum(log(eigen(param_value(s, eta), only.values = TRUE)$values)))
#> [1] TRUE

# A rank-deficient family declares the deficiency at construction, and the
# null basis really is annihilated by the matrix at any free value.
r <- scaled_matrix(diag(c(1, 1, 0, 0)))
c(rank = r@rank, null_columns = ncol(r@null_basis))
#>         rank null_columns 
#>            2            2 
max(abs(param_value(r, 0.6) %*% r@null_basis))
#> [1] 0
```
