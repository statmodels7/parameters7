# The Abstract Class of a Constrained Parameter

A parameter is a map from an unconstrained vector \\\eta \in
\mathbb{R}^{d}\\ onto a quantity that lives in a constrained set: a
symmetric positive definite matrix, a probability vector, a stochastic
matrix. `parameter` is the abstract S7 class that map belongs to. It
records how many free values the map consumes, what each of them is
called, and whatever the family needs in order to evaluate itself.
Because the domain is the whole of \\\mathbb{R}^{d}\\, an optimizer may
move \\\eta\\ anywhere and every point it visits still maps to a valid
value, so no constraint has to be policed during a fit.

The class is abstract and is not meant to be instantiated. Call a
constructor:
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
for an unstructured covariance,
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md) or
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
for a structured one,
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
for a probability vector,
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
for a Markov chain.

## Usage

``` r
parameter(
  param_name = character(0),
  n_free = integer(0),
  free_names = character(0),
  param_params = list()
)
```

## Arguments

- param_name:

  A single character string naming the family, used in the error
  messages the validators raise and in the object's `print` output.
  `"log_cholesky"`, `"ar1"` and so on.

- n_free:

  The length \\d\\ of the free vector: a single non-negative integer.
  Zero is legal and describes a family with nothing to estimate. The
  validator rejects a vector, a negative value, or a length that
  disagrees with `free_names`.

- free_names:

  A character vector of length `n_free`, one label per free value, in
  the order the free vector holds them. Fixed at construction and part
  of the interface: consumers build their parameter tables from these
  labels, so the ordering is not free to change. The validator rejects a
  duplicated label and a length other than `n_free`.

- param_params:

  A list of whatever the family needs in order to evaluate itself, read
  only by that family's own methods.
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  stores the row and column index of each free value here;
  [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
  stores its component matrices. Nothing outside the family looks inside
  it.

## Value

An object of class `parameter`, with properties

- `param_name`:

  character, the family name.

- `n_free`:

  integer, the length \\d\\ of the free vector.

- `free_names`:

  character of length `n_free`.

- `param_params`:

  list, the family's own data.

The class is abstract, so a useful object comes from a constructor such
as
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
and carries that constructor's subclass. A matrix family returns a
[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
which adds `dimension`, `rank` and `null_basis`.

## A parameter owns its dimension

`log_cholesky(3)` and `log_cholesky(4)` are different objects, holding
\\d = 6\\ and \\d = 10\\ free values. Fixing the size at construction
means `n_free` and `free_names` can be answered before any data exist,
and both are needed then: a model builds its parameter table from the
names while it is still assembling the design.

## What a subclass must supply

[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
is the one compulsory method. Given it, a subclass inherits
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md),
each a product stencil applied to the map itself. A closed form
registered later replaces the inherited numerical method through
dispatch, and no caller changes.

Four generics are not inherited here:

- [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md),
  the inverse map, has a base method that throws and names the family.
  An inverse found by optimization would hand back a plausible \\\eta\\
  for a value that is outside the set the family parametrizes, so the
  inverse is either written out exactly or refused.

- [`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md),
  [`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
  and
  [`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
  are registered on
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  not on this class. A probability vector has no log-determinant, and
  asking for one fails at dispatch instead of returning a number.

## Free names label the coordinate, not the quantity

The families here follow one convention. Where a link carries a
constrained quantity onto the free scale, the label records that link: a
variance appears as `"log_scale"` and a correlation as `"z_rho"`. Where
the coordinate is already unrestricted the label is the plain name of
the quantity, as the below-diagonal entries `"L2.1"` of a Cholesky
factor are.

The distinction matters outside the family. A consumer flattens the free
vector into scalar parameters carrying identity links, so a label
promising a bounded quantity would report a number that is not on that
scale: a free value of 1.4 called `rho` reads as a correlation of 1.4,
while called `z_rho` it reads as the correlation \\\tanh(1.4) = 0.885\\.

## Notation

\\\eta\\ is the **free vector**, the point on the unconstrained scale
that an optimizer moves, and \\d\\ its length, the object's `n_free`.
Note that \\\eta\\ means a linear predictor elsewhere in the toolkit;
here it is always the free vector of a parametrization. \\p\\ is the
side of the matrix a matrix family produces.

## See also

[`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
for the symmetric matrix branch.
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
and
[`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
for the map and its inverse,
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
for the derivative arrays,
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
for the quantities a family reports, and
[`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
to verify a parametrization written from scratch.

The constructors:
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md),
[`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md),
[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
[`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
[`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md),
[`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
[`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
[`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
[`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md),
[`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md),
and the compositions
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md),
[`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
[`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md).

## Examples

``` r
# Every parametrization in the package is a `parameter`.
s <- log_cholesky(3)
S7::S7_inherits(s, parameter)
#> [1] TRUE

# The four properties the class itself holds.
s@param_name
#> [1] "log_cholesky"
s@n_free
#> [1] 6
s@free_names
#> [1] "log_L1" "log_L2" "log_L3" "L2.1"   "L3.1"   "L3.2"  

# The dimension is fixed at construction, so the length of the free vector
# is known before any data arrive: p(p+1)/2 for an unstructured covariance.
d <- vapply(2:5, function(p) log_cholesky(p)@n_free, integer(1))
rbind(p = 2:5, n_free = d, `p(p+1)/2` = (2:5) * (3:6) / 2)
#>          [,1] [,2] [,3] [,4]
#> p           2    3    4    5
#> n_free      3    6   10   15
#> p(p+1)/2    3    6   10   15

# A family whose value is not a symmetric matrix inherits `parameter`
# directly, so it has no log-determinant to be asked for.
q <- simplex(3)
c(parameter = S7::S7_inherits(q, parameter),
  matrix_parameter = S7::S7_inherits(q, matrix_parameter))
#>        parameter matrix_parameter 
#>             TRUE            FALSE 
sum(param_value(q, c(0.7, -0.3)))          # a probability vector: sums to 1
#> [1] 1

# The free names record the link, so a free value of 1.4 in a `z_rho`
# coordinate is a correlation of tanh(1.4), not a correlation of 1.4.
ar1(4)@free_names
#> [1] "log_scale" "z_rho"    
param_value(ar1(4), c(0, 1.4))[1, 2]
#> [1] 0.8853516
tanh(1.4)
#> [1] 0.8853516
```
