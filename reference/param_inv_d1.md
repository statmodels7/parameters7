# Derivatives of the Inverse of a Matrix Parameter

`param_inv_d1()` returns the first derivatives of \\M^{-1}\\ in the free
values, and `param_inv_d2()` the second, where \\M\\ is the matrix that
[`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
returns.

## Usage

``` r
param_inv_d1(s, eta, ...)

param_inv_d2(s, eta, ...)
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
names its own. `param_inv_d2()`: a list with one matrix per pair in
[`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
at order 2, named as
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
names its own.

## Details

The default method differentiates the inverse through the derivatives of
\\M\\ itself, \$\$\partial_k M^{-1} = -M^{-1} A_k M^{-1}, \qquad
\partial\_{kl} M^{-1} = M^{-1}(A_l M^{-1} A_k + A_k M^{-1} A_l -
A\_{kl}) M^{-1},\$\$ with \\A_k\\ and \\A\_{kl}\\ from
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md).

Where \\M\\ is nearly singular along a direction that is not a
coordinate axis, \\M^{-1}\\ has large entries of nearly rank one, the
products above cancel, and the default method loses accuracy as the
condition number of \\M\\ grows.

The method for
[`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
is exact there. With \\M = L L^\top\\, \\G = L^{-1}\\, \\B_k =
G\\\partial_k L\\ and \\C_k = B_k + B_k^\top\\, \$\$\partial_k M^{-1} =
-G^\top C_k G, \qquad \partial\_{kl} M^{-1} = G^\top\big(B_l^\top C_k +
C_k B_l - \partial_l C_k\big) G,\$\$ where \\\partial_l B_k = -B_l B_k +
G\\\partial\_{kl} L\\. These products involve no cancellation: the
entries of \\G\\ grow as \\1/L\_{ii}\\ and each term is of the size of
the result.

## See also

[`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
for the inverse itself,
[`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
and
[`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
for the derivatives of \\M\\.

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
```
