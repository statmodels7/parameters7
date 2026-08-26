# The Scale, the Partial Autocorrelations and the Coefficients

Declares three groups: the marginal variance, the \\q\\ partial
autocorrelations that parametrize the family, and the \\q\\
autoregressive coefficients they produce. It is the only family in the
package that declares more quantities than it has free values, and the
reason is the third group.

## Arguments

- s:

  An
  [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  object.

- eta:

  A numeric vector of free values, of length \\q+1\\.

- ...:

  Ignored.

## Value

A list as described in
[`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md).

## Why the coefficients are declared

\\\phi_1, \dots, \phi_q\\ are what an autoregression is reported in, and
they appear nowhere in the covariance the fit prints, nor among the free
values, which are partial autocorrelations. Reporting a coordinate in
their place gives a reader the wrong number under the right name: at \\q
= 2\\ with free values \\(0, 0.9, -0.4)\\ the coefficients are \\(0.988,
-0.380)\\ while the partial autocorrelations are \\(0.716, -0.380)\\, so
\\\phi_1\\ and \\\rho_1\\ differ in the first decimal. At \\q = 1\\ they
coincide exactly, the first partial autocorrelation being the AR(1)
coefficient, which is why the distinction is invisible until \\q\\
passes 1.

## The Jacobian is not diagonal

A coefficient reads the **whole** chart, so its row of the Jacobian is
dense in the partial autocorrelation columns. Those derivatives come out
of the Levinson-Durbin recursion the family already propagates
derivative arrays through, so they are read off the first-order block
instead of recomputed. The variance's row is a single entry, and each
partial autocorrelation's is the link's own derivative.

## The coefficients' intervals are on the identity scale

The stationary region in the coefficients is not a box: it is bounded by
the roots of \\1 - \phi_1 z - \cdots - \phi_q z^q\\ lying outside the
unit circle, which no scalar transformation of a single coefficient
expresses. An interval that respected the constraint would not be an
interval, so the identity scale is the honest choice and a reported
coefficient interval may extend outside the stationary region. The
partial autocorrelations, whose chart **is** a box, get `"atanh"` and
stay inside it.
