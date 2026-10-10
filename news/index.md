# Changelog

## parameters7 0.23.0

- [`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
  of a
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  family takes its first and second derivatives from
  [`param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
  and
  [`param_inv_d2()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md),
  which are exact where the ordered-block-partition sum cancels. At
  `log L22 = -11.5` of a two-dimensional chart the sum gave the second
  derivative with a relative error of 4112 (wrong sign) and the first
  with 7.2e-7; both now agree with
  [`param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
  and
  [`param_inv_d2()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
  exactly. Orders three and four, and every other inner family, still
  use the partition sum.

- [`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
  and
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
  evaluate the log-determinant, its four derivative orders and the solve
  from the free value of the correlation, not from the correlation
  itself, which rounds to one of its bounds at large free values. At p =
  4 the log-determinant of
  [`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
  was wrong by 1.9e-7 at `logit_rho = 20` and stayed at -106.74 from 37
  on (the true value is -117.75 at 40); the log-determinant of
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
  stayed at -106.05 from `z_rho = 20` on (the true value is -235.84 at
  40). Both now agree with 60-digit references to 1e-12 relative at free
  values up to 40, at every order and in the inverse. The internal
  [`cs_logdet_chain()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_chain.md)
  and
  [`ar1_logdet_chain()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_chain.md)
  hold the closed forms; `cs_logdet_terms()` and `ar1_logdet_terms()`
  are removed.

- [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
  [`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
  and
  [`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
  evaluate every factor `1 - r^2` from the free value of the partial
  autocorrelation: the log-determinant of
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
  its derivatives and the innovation variances of its solve through
  [`sech2()`](https://statmodels7.github.io/parameters7/reference/sech2.md)
  and
  [`log_sech2()`](https://statmodels7.github.io/parameters7/reference/sech2.md),
  and the factor `1 / (1 - r^2)` of the two precisions as `cosh(z)^2`
  through
  [`cosh2_derivs()`](https://statmodels7.github.io/parameters7/reference/cosh2_derivs.md).
  Formed from `r`, they lost their digits as `r` rounded to one of its
  bounds: from `z = 20` on, the derivatives of
  [`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
  and
  [`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
  had a relative error of 1, and the log-determinant of
  `autoregressive(4, 1)` stayed at -106.05. The derivative arrays of
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
  and
  [`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
  also need linkfunctions7 0.5.2, whose rhobit inverse derivatives keep
  their accuracy for a large free value; it is now the minimum version.
  All of these agree with 90-digit references to 1e-12 at free values up
  to 20, and those of
  [`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
  up to 40. `w_derivs()` and `log_affine_derivs()` are removed.

- The default third and fourth log-determinant derivatives accept a
  family that writes its own
  [`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
  which they difference, whatever its
  [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
  and
  [`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md);
  they rejected it.

- [`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md)
  symmetrizes its estimates for a matrix family, as the higher orders
  do.

- The derivatives of
  [`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
  have a diagonal of exact zeros at every order; it was zero only up to
  rounding.

- [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md):
  the log-determinant derivatives of a rank-deficient family use the
  pseudo-inverse on the complement of the declared null space; they
  stopped with a LAPACK error. The constructor rejects linearly
  dependent components, whose weights would not be identified, and
  [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
  names its result by `free_names`.

- [`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md),
  which decides positive definiteness for the inverse maps of
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  and
  [`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
  and for the default factor, uses a relative eigenvalue tolerance of
  1e-14 instead of 1e-12. Over 20000 round trips of
  `correlation_matrix(4)` at free values of standard deviation 2.5, the
  positive definite matrices rejected fall from 54 to 5.

- [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  compares the first and second derivatives and the log-determinant
  Hessian at every free vector of its sweep, not only at the three
  constant ones; an asymmetric value fails the membership check; and for
  a family that is not a matrix, a derivative order from a numerical
  fallback is reported as NOT CHECKED, and the derivative thresholds
  scale with `tol`.

- [`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md)
  accepts a block without free values, and its
  [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
  names its result.
  [`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)’s
  [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
  names its result by `free_names`. `scaled_matrix(link = NULL)` returns
  empty results at derivative orders three and four, as at one and two.
  The inverse maps of
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
  [`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
  and
  [`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
  reject a value outside the range of the link.
  [`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md)
  rejects a non-numeric, logical, infinite or too large `m`, the inverse
  maps reject a matrix with an infinite entry, and the class validators
  reject a non-scalar `n_free`, `dimension` or `rank` with their own
  messages.

- The documentation is rewritten in a plain register, and its statements
  are checked against the code: the counts of families and wrappers, the
  claims on accuracy and cost, and the descriptions of the base-class
  methods.

## parameters7 0.22.1

- The test of
  [`param_inv_d2()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
  at a nearly singular covariance compares it with symbolic second
  derivatives of the inverse, to 1e-12. The Richardson difference that
  it used had an error of about 1e-6 there, above the tolerance on Linux
  (1.05e-6); the method agrees with the symbolic values to 4e-16.

## parameters7 0.22.0

- Every helper forms only the orders that its caller reads. The value of
  a parameter, its solve and its log-determinant no longer evaluate any
  derivative, and a derivative of order k no longer forms the orders
  above k. `compose4()`, which always composed four orders, is replaced
  by
  [`compose_order()`](https://statmodels7.github.io/parameters7/reference/compose_order.md),
  one order per call;
  [`power_derivs()`](https://statmodels7.github.io/parameters7/reference/power_derivs.md),
  `log_affine_derivs()`,
  [`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md),
  [`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md),
  [`dr_scale_derivs()`](https://statmodels7.github.io/parameters7/reference/dr_scale_derivs.md)
  and
  [`sum_struct_weight_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_weight_derivs.md)
  take the highest order wanted, and
  [`linkinv_upto()`](https://statmodels7.github.io/parameters7/reference/linkinv_upto.md)
  replaces the five-element link tables that each family built. The
  compiled Levinson-Durbin recursion of
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
  propagates only the tensors up to the requested order, and the value
  is the recursion on numbers alone. The derivatives are unchanged, and
  the record of order k is the leading columns of the order-four record.
  At p = 40, the value of an `autoregressive(order = 3)` takes 0.20 ms
  against 3.55 and its first derivative 0.8 against 4.0; the value of an
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
  1.3 ms against 6.3.

- [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
  has one compiled kernel per derivative order (`ar_d1_cpp()` to
  `ar_d4_cpp()`, and `ar_value_cpp()`), which differentiates the
  Levinson-Durbin recursion in the partial autocorrelations, keeps one
  entry per multiset of indices, applies the links once through the
  partial Bell polynomials and returns the components of that order
  only. The Toeplitz matrices are filled in compiled code
  (`ar_toeplitz_cpp()`). At q = 4, p = 200:
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  0.15 ms against 0.60,
  [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
  0.43 against 2.8,
  [`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
  15.5 against 55.6. The former recursion stays compiled as
  `ar_taylor_jet_cpp()`, the independent reference of the tests, which
  the new kernels match to 1e-12 at every order.

## parameters7 0.21.0

- [`param_inv_d1()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md)
  and
  [`param_inv_d2()`](https://statmodels7.github.io/parameters7/reference/param_inv_d1.md),
  new generics: the first and second derivatives of the inverse of a
  matrix parameter in its free values. The default method, for any
  `matrix_parameter`, is `-M^-1 A_k M^-1` and its second-order
  counterpart. The method for
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  avoids the cancellation of that expression: with `G = L^-1`,
  `B_k = G dL_k` and `C_k = B_k + B_k'`, the first derivative is
  `-G' C_k G` and the second `G'(B_l' C_k + C_k B_l - dC_k/dl)G`.

  Where the matrix is nearly singular along a direction that is not a
  coordinate axis, the default method loses accuracy: at the free value
  `log L22 = -11.5` of a two-dimensional log-Cholesky chart it gives the
  first derivative with a relative error of 5.6e-8 and the second one
  entirely wrong, where the method for
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  agrees with a Richardson difference of the inverse to 5e-12. At
  ordinary values the two agree to 1e-15. The multivariate families of
  distributions7 read these for the derivatives of `Sigma^-1`.

## parameters7 0.20.0

- [`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
  writes its derivative arrays out instead of taking them from the
  ordered-block-partition sum of
  [`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md).
  The value comes from the prediction form `Omega = U' diag(tau) U`. The
  lower-order rows of U are the coefficients of the same family at that
  order, so their derivative arrays come from the compiled
  Levinson-Durbin recursion run once per order, and a component
  differentiating in a partial autocorrelation that an order does not
  reach is exactly zero. `tau` is a product of one univariate factor per
  free value, so a mixed derivative of it is a product of univariate
  derivatives. The rest is the Leibniz rule over three factors, taken
  twice so that a component costs 2^m matrix products and not 3^m. At
  order four against the composition: 11.2 times faster at p = 6, q = 2,
  11.4 at p = 20, q = 3 and 22.8 at p = 100, with the two agreeing to
  7e-11 over six shapes and two free vectors each.

- Sub-multisets are held as count vectors indexed by an integer code, in
  place of string keys built with
  [`sort()`](https://rdrr.io/r/base/sort.html) and
  [`paste()`](https://rdrr.io/r/base/paste.html); the subset structure
  of an order is memoized on the object; and the lower-order families
  are built once by the constructor. A fourth-order component at p = 6,
  q = 2 goes from 19.21 ms to 3.74 ms, every number unchanged.

- `w_derivs()` gives the derivatives of `1/(1 - rho^2)`, exact at every
  order by partial fractions.

## parameters7 0.19.0

- The `role` property and its argument are removed. A parametrization no
  longer records the side of a model on which it is used; the consumer
  states that, and `distributions7`’s `mvgaussian1_distrib()` and
  `mvgaussian2_distrib()` state it by which family they are.
  `distributions7` never read the label, `modelterms7::random()` reads
  the side from the distribution’s `params_interpretation`, and its one
  use, a guard in `penalties7::structured_penalty()`, is now written
  with
  [`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md).

- [`inverse_of()`](https://statmodels7.github.io/parameters7/reference/inverse_of.md)
  is the fifth composition wrapper: the family whose value is another
  family’s inverse, carrying the inner family’s free vector unchanged.
  Its derivatives are the sum over the **ordered** set partitions of the
  differentiated positions,
  `d^I N = sum over (B_1, ..., B_q) of (-1)^q N (d^B_1 S) N ... (d^B_q S) N`,
  the blocks being ordered because matrices do not commute; the number
  of terms is the Fubini number of the order, 1, 3, 13 and 75. The
  log-determinant is the inner one negated, and
  [`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
  is the inner family’s value. It passes
  [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  on the eight families on which it was tried, including the three that
  are not closed under inversion. A rank-deficient family is rejected.

- [`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
  and
  [`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
  are the two families whose inverse has a structure of its own: the
  tridiagonal precision of a first-order autoregression and the banded
  precision of an order-q one. Both carry the free vector of the family
  that they invert.
  [`ar1_inv()`](https://statmodels7.github.io/parameters7/reference/ar1_inv.md)
  writes its derivative arrays out: the value is a function of the scale
  times a matrix function of the correlation, and the three distinct
  entries are written in `w = 1/(1 - rho^2)`, whose derivatives are
  exact by partial fractions. At p = 20 and order four it takes 0.44 ms
  against 28.66 ms for the composition wrapper, and the two agree to
  8e-14.
  [`autoregressive_inv()`](https://statmodels7.github.io/parameters7/reference/autoregressive_inv.md)
  takes the closed value, log-determinant and solve of the family that
  it inverts, and its derivative arrays from the composition.

## parameters7 0.18.0

- [`check_positive_link()`](https://statmodels7.github.io/parameters7/reference/check_positive_link.md)
  tests both ends of the link, and every family taking one rejects a
  link defined on part of the real line. The check read only the lower
  bound of the link’s range, so `sqrt_link()`, `inverse_link()`,
  `inverse_sq_link()` and `power_link()` at a positive exponent were
  accepted, although they reach the positive entries from the positive
  predictors alone. Under `sqrt_link()` the map is not injective: the
  first coordinate of
  `param_free(autoregressive(5, 1, link_scale = sqrt_link()), M)` came
  back with its sign flipped, and
  [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  stopped with an error that did not name the link. The second condition
  is read through
  [`linkfunctions7::eta_bounds()`](https://statmodels7.github.io/linkfunctions7/reference/eta_bounds.html).
  The eight constructors taking a link are covered:
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md),
  [`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md),
  [`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md),
  [`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md),
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md),
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md),
  [`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
  and
  [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md).

## parameters7 0.17.0

- [`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
  rejects a `correlation` block that does not have a unit diagonal. A
  block carrying a scale of its own was accepted, and the composite then
  had one free value too many: the scale moves between D and R without
  changing the matrix, so
  [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
  returned a different point from the one that
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  was given. The diagonal is read at two probe free vectors,
  `0.3, 0.4, ...` and a constant `-0.4`, fixed and not drawn so that a
  rejection is reproducible; the zero vector is not used, because a
  scale-carrying family on a log link has a unit diagonal there.

## parameters7 0.16.0

- [`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
  and
  [`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md)
  signal an error where the family supplies no analytic
  [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
  or
  [`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md).
  Both fallbacks difference
  [`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md),
  which is an exact identity given those two arrays, so the differencing
  is one numerical layer where they are written out and two where they
  are not. On a 4 by 4 AR(1) covariance whose family supplies
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  alone, order three came back 7.5e-3 against a quantity of size 4.45
  and order four 9.07 against a quantity of size 2.17. The families of
  the package are not affected, and the message names the missing
  method.

## parameters7 0.15.0

- `scaled_dlog()` is removed, and the two log-determinant methods that
  called it read
  [`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md),
  which computes the same quantity at every order from one to four. The
  removed function answered order 2 with the fourth-order expression (on
  a square-root link at a free value of 1.5 it returned -2.37037 where
  the second derivative of the log inverse link is -0.888889); its two
  callers passed only 3 or 4, so no reported number changes.

## parameters7 0.14.0

- The four composition wrappers then available
  ([`kron_identity()`](https://statmodels7.github.io/parameters7/reference/kron_identity.md),
  [`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md),
  [`dr_prod()`](https://statmodels7.github.io/parameters7/reference/dr_prod.md)
  and
  [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md))
  label the rows and columns of the value and of the four derivative
  orders, the convention that
  [`name_dims()`](https://statmodels7.github.io/parameters7/reference/name_dims.md)
  states; they returned matrices without dimnames.

- The shapes-and-names check of
  [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  reads the margins of
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md),
  [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
  and
  [`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
  as well as the declared names.

- [`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
  and
  [`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
  on a
  [`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
  keep returning matrices without dimnames: the convention covers the
  value and the four derivative orders.

## parameters7 0.13.0

- [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  leaves the caller’s random stream as it was. The matrix battery drew
  its random free vectors from the caller’s stream, and the branch for a
  family that is not a matrix called
  [`set.seed()`](https://rdrr.io/r/base/Random.html) and replaced the
  caller’s state (`set.seed(42); rnorm(1)` gave 1.370958, and -1.172560
  after `check_parameter(simplex(3))`). Both branches now draw from a
  fixed seed and restore the `.Random.seed` that they found on entry
  through [`on.exit()`](https://rdrr.io/r/base/on.exit.html). The
  internal helpers are
  [`capture_seed()`](https://statmodels7.github.io/parameters7/reference/capture_seed.md)
  and
  [`restore_seed()`](https://statmodels7.github.io/parameters7/reference/capture_seed.md).

## parameters7 0.12.0

- `correlation_matrix(1)` builds a constant. A one by one correlation
  matrix has no angles, so `n_free` is 0; the constructor reported 1,
  with the free name `"z."`, and
  [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md),
  [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
  and
  [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  stopped with an error.
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  now returns the identity at `numeric(0)`, the derivative lists are
  empty, the log-determinant is 0, and
  [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  passes with its two log-determinant derivative rows reported as NOT
  CHECKED. The family can therefore be composed inside
  [`block_diag()`](https://statmodels7.github.io/parameters7/reference/block_diag.md).

## parameters7 0.11.0

- The numerical fallbacks take their stencils from numericals7:
  [`mixed_stencil()`](https://statmodels7.github.io/parameters7/reference/mixed_stencil.md)
  reads the offsets and weights from `fd_weights()` in place of a
  transcribed table,
  [`fd_step()`](https://statmodels7.github.io/parameters7/reference/fd_step.md)
  forwards to numericals7, and
  [`fd_along()`](https://statmodels7.github.io/parameters7/reference/fd_along.md)
  assembles every difference along a free value. Every number is
  unchanged.

## parameters7 0.10.0

- `sum_struct(components, link)`: a non-negative combination of fixed
  symmetric positive semidefinite matrices, which is the
  variance-components covariance and also the matrix that a penalty with
  one smoothing parameter per component assembles. The value is linear
  in the weights, so a derivative component is zero unless every index
  names the same free value. The log-determinant is not separable, and
  its derivatives come from the cyclic trace expansion
  `(-1)^(n-1) sum tr(M^-1 P_s1 ... M^-1 P_sn)`, the sum running over the
  (n-1)! orderings of the last n-1 indices, counted with multiplicity,
  and carried onto the free scale by a chain rule with a diagonal
  Jacobian. Counting the coinciding orderings once would make a third
  derivative in one weight too small by a factor of 2 and a fourth by 6.

  The rank is fixed at construction from the components stacked and
  individually normalized, the null space of a sum of positive
  semidefinite matrices being the intersection of theirs. On second and
  first differences over five coefficients, the null residual stays at
  1e-16 with the weights fourteen orders of magnitude apart.

## parameters7 0.9.0

- `block_diag(...)`: a block diagonal of distinct matrix parameters,
  each carrying its own stretch of the free vector. Every component
  whose indices do not all belong to one block is exactly zero, at every
  order and for the log-determinant as well as for the value; the rank
  and the null basis come from the blocks.

- `dr_prod(dimension, correlation, link)`: a covariance written as
  `D R D`, whose coordinates are the linked standard deviations followed
  by those of the correlation. Every derivative factorizes as
  `d^S_D (d_i d_j) * d^S_R R_ij`, the two groups of free values being
  disjoint, and the log-determinant `2 sum log d_j + log|R|` is
  separable in the scales and from the correlation. The correlation
  block must have full rank, because a deficient one would give the
  product the null space `D^-1 ker R`, which moves with the free vector.

## parameters7 0.8.0

- `kron_identity(structure, m)`: m identical diagonal blocks of one
  matrix parameter sharing a free vector, every quantity a linear lift
  of the inner one. The first composition wrapper, for grouped random
  effects.

## parameters7 0.7.0

- The log-Cholesky derivative assembly is compiled. Every derivative of
  the factor is a single-entry matrix, so each Leibniz term is one row,
  one column or one cell; `chol_leibniz_cpp` uses that and serves all
  four orders. At p = 8: order 4 from 8.44 s to 0.28 s (the 82,251
  result matrices being most of what remains), order 2 from 6 ms to 1.2
  ms. The dense R assembly stays as the twin `.chol_leibniz_r`, compared
  with the compiled one in the tests.

## parameters7 0.6.0

- [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
  no longer consumes jets: the Levinson-Durbin recursion propagates its
  derivative arrays in compiled code (`ar_taylor_cpp`), each tracked
  quantity holding its value and its symmetric tensors to fourth order,
  combined by the product rule written out per order. Rcpp enters
  Imports and LinkingTo.
  [`ar_prediction()`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md)
  and the log-determinant read the link inverses directly. The q = 1
  derivatives are tested against products of link derivatives, formulas
  that share no code with the kernel.

## parameters7 0.5.0

- The jets move to numericals7.
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
  consumes them through numericals7, `jet_layout()`, `jet_var()` and
  `set_partitions()` are no longer exported from here, and
  [`tuple_indices()`](https://statmodels7.github.io/parameters7/reference/tuple_indices.md)
  forwards to the enumeration in numericals7.

## parameters7 0.4.0

- `jet_compose()` applies a smooth scalar function to a jet given that
  function’s own five derivatives, summing over the set partitions of
  the positions of an index tuple, so `exp`, `log`, an arbitrary power,
  `sqrt`, `gamma`, `lgamma`, `digamma` and `trigamma` each need five
  numbers. `Ops` and `Math` dispatch on the class, so a map written as
  `mu / gamma(1 + 1 / sigma)` carries every partial derivative to fourth
  order. Comparison operators and the non-smooth functions are rejected.

- `jet_layout()`, `jet_var()` and `set_partitions()` are exported, for
  seeding jets and enumerating partitions. distributions7 uses them in
  `reparametrize()`.

## parameters7 0.3.0

- [`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
  declares the quantities that a family describes, with the Jacobian of
  the map from the free vector and the scale on which each interval is
  built, so that a consumer can report them by the delta method.
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
  and
  [`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
  declare a marginal variance and a correlation,
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
  adds the coefficients produced by the Levinson-Durbin recursion,
  [`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
  its multiplier and
  [`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
  its probabilities. The base class declares nothing.

- Free names name the coordinate and not the quantity that it produces.
  Where a link carries a constrained quantity onto the free scale the
  name records that link:
  [`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
  reports `log_scale`,
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
  `log_scale` and `z_rho`,
  [`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
  `logit_rho`,
  [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
  `z_pacf1` and so on, and
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  `log_d1` under its default link. Names built on a coordinate that is
  already unrestricted are unchanged: `log_L1` and `L2.1`, `S2.1`,
  `alr1`, `z2.1`. The previous names suggested bounded quantities while
  reporting values on the free scale.

- `autoregressive(dimension, order)` carries the covariance of a
  stationary autoregression of any order, parametrized by the partial
  autocorrelations, each in `(-1, 1)` and carried onto the coefficients
  by the Levinson-Durbin recursion (Barndorff-Nielsen and Schou, 1973).
  The log-determinant is `p log(g0) + sum (p-k) log(1 - r_k^2)` and the
  precision is banded with the order as bandwidth.
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
  stays as the order-one case written out.

- [`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
  carries a correlation matrix in the spherical parametrization of
  Rapisarda, Brigo and Mercurio (2007), and its log-determinant is twice
  the sum of the logarithms of the sines.

- [`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
  and
  [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md):
  two free values whatever the dimension, a scale and a correlation,
  with a closed log-determinant separable in the two free values and a
  closed inverse (Sherman-Morrison for compound symmetry, tridiagonal
  for AR(1)). Compound symmetry bounds its correlation below at
  `-1/(p-1)`, where the matrix stops being definite.

- All three are closed form to fourth order, in the value and in the
  log-determinant, through `compose4()` for the chain rule and
  [`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md)
  for the derivative of a Gram product.

## parameters7 0.2.0

- Renamed from covstructs7, with the API renamed to match: `struct_*`
  became `param_*`, the base class is `parameter`, and the symmetric
  positive semidefinite branch is the class `matrix_parameter`, which
  holds the rank, the log-determinant, the solve and the factor.

- Derivatives to fourth order:
  [`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md),
  [`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md),
  [`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
  and
  [`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.md),
  closed form for every family of the package and served by
  single-stencil numerical fallbacks otherwise.

- [`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
  carries a probability vector in the additive log-ratio chart;
  [`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
  is its row-wise extension to row-stochastic matrices;
  [`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
  is the matrix logarithm chart on the positive definite cone, with a
  linear log-determinant, an exact inverse, and Frechet derivatives by
  Daleckii-Krein with divided differences computed by the Opitz theorem.

## parameters7 0.1.0

- First release. `parameter` is the base class of constrained matrix
  parameters: a map from an unconstrained vector to a symmetric matrix,
  with its first and second derivatives, its log-determinant and its
  solves. Only
  [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  is compulsory; the rest has methods on the base class.

- Four families.
  [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  is the unstructured positive definite case in the parametrization of
  Pinheiro and Bates (1996), whose log-determinant is linear in the free
  vector.
  [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  and
  [`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
  carry their entries through a `linkfunctions7` link.
  [`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
  is a fixed matrix carried by one scale, covering the ridge, the spline
  penalties with one smoothing parameter and the fully known case.

- Rank-deficient matrices are accepted, as a spline penalty requires.
  [`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
  returns the log pseudo-determinant there, and its derivative under a
  scaled parameter is the rank.

- The rank and the null space are computed once at construction, from
  the components, by
  [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md).
  On a tensor-product penalty of rank 28 out of 32, counting the
  eigenvalues of the assembled sum above a relative tolerance gives 28
  at equal scales and 16 at a ratio of 1e12.

- [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  runs nine checks, each against a route that the implementation does
  not take, and reports a quantity from the base class as not checked.
