# Construct a Scalar Multiple of the Identity

Returns an object holding the map \\M = \tau I\\ with \\\tau =
h(\eta_1)\\ positive: a diagonal matrix with **one** free value shared
by every entry, whatever \\p\\ is. It is the simplest parametrization in
the package, and it is what a random effect with a single variance
component needs.

## Usage

``` r
scalar_matrix(
  dimension,
  link = linkfunctions7::log_link(),
  role = c("either", "covariance", "precision")
)
```

## Arguments

- dimension:

  The side \\p\\ of the matrix. A single positive whole number, finite
  and at least 1; anything else throws.

- link:

  A linkfunctions7 link carrying the free value onto the positive scale,
  [`linkfunctions7::log_link()`](https://statmodels7.github.io/linkfunctions7/reference/log_link.html)
  by default. It must map onto the positive half line, so
  `identity_link()` is rejected, and from the whole real line, which
  rules out `sqrt_link()` and its relatives; see
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  for the two conditions.

- role:

  A label recording which side of a model the matrix parametrizes:
  `"either"` (the default), `"covariance"` or `"precision"`. No numeric
  result depends on it.

## Value

An object of class
[`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md),
with `n_free` 1, `free_names` a single link-tagged label, `rank` equal
to `dimension`, an empty `null_basis`, `param_name` `"scalar"`, and
`param_params` holding `link` and `shared = TRUE`.

## One free value, whatever the dimension

`n_free` is 1 at every `dimension`, and `free_names` is `log_scale`
under the default link. The multiplicity shows up in the
log-determinant, which is \\p \log h(\eta_1)\\, so the gradient carries
a factor of \\p\\ and a step of one in the free value moves the
log-determinant by \\p\\ under the log link.

## The same matrix from the other direction

Under the default log link this is `scaled_matrix(diag(dimension))`, a
positive multiple of a fixed matrix that happens to be the identity. The
two give the same value at the same free value. Prefer this one when the
matrix is the identity and
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
when it is not, and note that
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
admits a rank-deficient fixed matrix while this family is always of full
rank.

## Notation

\\\eta_1\\ is the single free value, \\\tau = h(\eta_1)\\ the shared
positive entry, \\p\\ the side of the matrix and \\h = g^{-1}\\ the
inverse link.

## See also

[`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
for one free value per entry,
[`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
for a multiple of an arbitrary fixed matrix, and
[`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
to replicate any parametrization over independent groups.

## Examples

``` r
# One free value at any dimension.
s <- scalar_matrix(3)
c(n_free = s@n_free, dimension = s@dimension)
#>    n_free dimension 
#>         1         3 
s@free_names
#> [1] "log_scale"
param_value(s, log(2))
#>    v1 v2 v3
#> v1  2  0  0
#> v2  0  2  0
#> v3  0  0  2

# It is scaled_matrix() on the identity, reached from the other side.
all.equal(param_value(s, log(2)),
          param_value(scaled_matrix(diag(3)), log(2)),
          check.attributes = FALSE)
#> [1] TRUE

# The shared value enters the log-determinant p times over, so the gradient
# carries a factor of p.
c(logdet = param_logdet(s, log(2)), p_log_2 = 3 * log(2))
#>   logdet  p_log_2 
#> 2.079442 2.079442 
param_dlogdet(s, log(2))
#> log_scale 
#>         3 
```
