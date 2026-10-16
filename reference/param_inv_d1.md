# Derivatives of the Inverse of a Matrix Parameter

`param_inv_d1()`, `param_inv_d2()`, `param_inv_d3()` and
`param_inv_d4()` return the derivatives of first to fourth order of
\\M^{-1}\\ in the free values, where \\M\\ is the matrix that
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns.

## Usage

``` r
param_inv_d1(s, eta, ...)

param_inv_d2(s, eta, ...)

param_inv_d3(s, eta, ...)

param_inv_d4(s, eta, ...)
```

## Arguments

- s:

  An object inheriting from class
  [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md),
  of full rank.

- eta:

  A numeric vector of length `s@n_free`, finite in every entry.

- ...:

  Passed to the method. The methods in this package do not read it.

## Value

`param_inv_d1()`: a list of `s@n_free` matrices, named as
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
names its own. `param_inv_d2()`, `param_inv_d3()` and `param_inv_d4()`:
a list with one matrix per tuple in
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
at that order, named as
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md),
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
name their own.

## Details

The default method differentiates the inverse through the derivatives of
\\M\\ itself, \$\$\partial_k M^{-1} = -M^{-1} A_k M^{-1}, \qquad
\partial\_{kl} M^{-1} = M^{-1}(A_l M^{-1} A_k + A_k M^{-1} A_l -
A\_{kl}) M^{-1},\$\$ with \\A_k\\ and \\A\_{kl}\\ from
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md).
At third and fourth order it uses the sum over the ordered set
partitions that
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
writes out, with the arrays of
[`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
and
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md).

Where \\M\\ is nearly singular along a direction that is not a
coordinate axis, \\M^{-1}\\ has large entries of nearly rank one, the
products above cancel, and the default method loses accuracy as the
condition number of \\M\\ grows, roughly as its square at the higher
orders.

The method for
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
is exact there. With \\M = L L^\top\\, \\G = L^{-1}\\, \\B_k =
G\\\partial_k L\\ and \\C_k = B_k + B_k^\top\\, \$\$\partial_k M^{-1} =
-G^\top C_k G, \qquad \partial\_{kl} M^{-1} = G^\top\big(B_l^\top C_k +
C_k B_l - \partial_l C_k\big) G,\$\$ where \\\partial_l B_k = -B_l B_k +
G\\\partial\_{kl} L\\. These products involve no cancellation: the
entries of \\G\\ grow as \\1/L\_{ii}\\ and each term is of the size of
the result.

At third and fourth order the same method differentiates \\L G = I\\ and
\\M^{-1} = G^\top G\\ by the Leibniz rule. For a multiset \\S\\ of free
values, \$\$\partial_S G = -G \sum\_{\emptyset \ne U \subseteq S}
(\partial_U L)(\partial\_{S \setminus U} G), \qquad \partial_S M^{-1} =
\sum\_{U \subseteq S} (\partial_U G)^\top (\partial\_{S \setminus U}
G),\$\$ the subsets being taken over the positions of \\S\\. Every
factor is triangular and every derivative of \\L\\ has a single non-zero
entry, so these products do not cancel either.

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
for the inverse itself,
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
to
[`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
for the derivatives of \\M\\, and
[`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
for the family whose value is the inverse.

## Examples

``` r
s <- log_cholesky(2)
eta <- c(-0.5, -11.5, 0.07)
d <- param_inv_d1(s, eta)
# against a difference of the inverse, which here is itself accurate
h <- 1e-6
num <- (param_solve(s, eta + c(0, 0, h)) - param_solve(s, eta - c(0, 0, h))) /
  (2 * h)
max(abs(d[[3]] - num)) / max(abs(num))
#> [1] 6.444149e-12
names(param_inv_d2(s, eta))
#> [1] "log_L1:log_L1" "log_L2:log_L2" "L2.1:L2.1"     "log_L1:log_L2"
#> [5] "log_L1:L2.1"   "log_L2:L2.1"  

# At third order the default route has lost the component that the exact
# one keeps: the inverse is quadratic in L2.1, so this one is zero.
sandwich <- S7::method(param_inv_d3, matrix_parameter)(s, eta)
c(exact = param_inv_d3(s, eta)[["L2.1:L2.1:L2.1"]][2, 2],
  sandwich = sandwich[["L2.1:L2.1:L2.1"]][2, 2])
#>        exact     sandwich 
#> 0.000000e+00 2.397528e+13 
```
