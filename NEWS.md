# parameters7 0.4.0

* Jets gain a composition rule and a vocabulary of functions, and can be
  written in ordinary R. `jet_compose()` applies a smooth scalar function to a
  jet given that function's own five derivatives, summing over the set
  partitions of the positions of an index tuple, so `exp`, `log`, an arbitrary
  power, `sqrt`, `gamma`, `lgamma`, `digamma` and `trigamma` each cost five
  numbers rather than a chain rule of their own. `Ops` and `Math` dispatch on
  the class, so a map reads `mu / gamma(1 + 1 / sigma)` and carries every
  partial derivative to fourth order with it. Comparison operators and the
  non-smooth functions are refused: a branch taken on a jet would keep one
  side's derivatives and report them as the whole expression's.

* `jet_layout()`, `jet_var()` and `set_partitions()` are exported, which is
  the seam a consumer needs to seed jets and to enumerate partitions.
  distributions7 uses it for `reparametrize()`.

# parameters7 0.3.0

* `param_readable()` declares the quantities a family is about, with the
  Jacobian of the map from the free vector and the scale each interval is
  built on, so that a consumer reports them by the delta method instead of
  reporting coordinates. `ar1()` and `compound_symmetry()` declare a marginal
  variance and a correlation, `autoregressive()` adds the coefficients that
  the Levinson-Durbin recursion produces -- read off the first-order component
  of jets it already computes -- `scaled_matrix()` its multiplier and
  `simplex()` its probabilities. The base class declares nothing. A
  multivariate fit in distributions7 prints the result as a block of its own.

* Free names name the coordinate rather than the quantity it produces, and
  every family now follows the one convention. Where a link carries a
  constrained quantity onto the free scale the name records that link, so
  `scalar_matrix()` reports `log_scale`, `ar1()` reports `log_scale` and
  `z_rho`, `compound_symmetry()` reports `logit_rho`, `autoregressive()`
  reports `z_pacf1` and so on, and `diagonal_matrix()` reports `log_d1`
  under its default link and `sqrt_d1` under a square-root one. Names built
  on a coordinate that is already unrestricted are unchanged: `log_L1` and
  `L2.1`, `S2.1`, `alr1`, `z2.1`. The previous names promised bounded
  quantities and were reported on the free scale, so a partial
  autocorrelation appeared as 0.97 with a confidence interval reaching past
  one.

* `autoregressive(dimension, order)` carries the covariance of a stationary
  autoregression of any order. The stationary region in the coefficients is
  not a box -- at order two it is already a triangle -- so the chart is the
  partial autocorrelations, each free in `(-1, 1)` and carried onto the
  coefficients by the Levinson-Durbin recursion, the transformation of
  Barndorff-Nielsen and Schou (1973). The map is polynomial, so its
  derivatives to fourth order come from carrying jets through the recursion
  rather than from expanding it; `jet_mul()` is the Leibniz rule again, in
  its scalar form. The log-determinant is `p log(g0) + sum (p-k) log(1-r_k^2)`
  and the precision is banded of the order's width. `ar1()` stays as the
  order-one case written out, and the two agree to machine precision.

* Three families for structured covariances. `correlation_matrix()` carries a
  correlation matrix in the spherical parametrisation of Rapisarda, Brigo and
  Mercurio (2007): the rows of the Cholesky factor are points on the unit
  sphere in angular coordinates, so the unit diagonal and the positive
  definiteness hold by construction, and the log-determinant is twice the sum
  of the logarithms of the sines.

* `compound_symmetry()` and `ar1()` are the economical families: two free
  values whatever the dimension, a scale and a correlation. Both have a closed
  log-determinant separable in the two free values, so every mixed derivative
  of it is exactly zero, and both return their inverse in closed form --
  compound symmetric again by Sherman-Morrison, tridiagonal for AR(1), the
  precision of a Markov process. Compound symmetry bounds its correlation
  below at `-1/(p-1)`, which is where the matrix stops being definite, rather
  than at `-1`.

* All three are closed form to fourth order, in the value and in the
  log-determinant. Two helpers carry that: `compose4()` for the chain rule to
  fourth order and `leibniz_gram()` for the derivative of a Gram product,
  which the log-Cholesky family now shares.

# parameters7 0.2.0

* Renamed from covstructs7, with the API renamed to match: the object is a
  constrained PARAMETER -- a map from an unconstrained vector onto the set
  where the parameter lives -- and not only a covariance. `struct_*` became
  `param_*`, the base class is `parameter`, and the symmetric positive
  semidefinite branch moved to the new abstract class `matrix_parameter`,
  which owns the rank, the log-determinant, the solve and the factor.

* Derivatives to FOURTH order. `param_d3()`/`param_d4()` and
  `param_d3logdet()`/`param_d4logdet()`, closed form for every shipped
  family and served by single-stencil numerical fallbacks otherwise.

* Three new families. `simplex()` carries a probability vector in the
  additive log-ratio chart, its derivatives closing over the value by the
  cumulant recursion of the categorical indicator; `transition_matrix()` is
  its row-wise extension to row-stochastic matrices; `matrix_log()` is the
  matrix logarithm chart on the positive definite cone, with the
  log-determinant linear, the inverse exact, and Frechet derivatives by
  Daleckii-Krein with divided differences computed by the Opitz theorem.

# parameters7 0.1.0

* First release. `parameter` is the abstract class of constrained matrix
  parameters: a map from an unconstrained vector to a symmetric matrix,
  together with its first and second derivatives, its log-determinant and the
  solves a likelihood asks of it. Only `param_value()` is compulsory;
  everything else has a method registered on the base class, so a new
  parameter is a subclass and one method.

* Four families. `log_cholesky()` is the unstructured positive definite case
  in the parametrisation of Pinheiro and Bates (1996), whose log-determinant
  is linear in the free vector. `diagonal_matrix()` and `scalar_matrix()` carry
  their entries through a `linkfunctions7` link, which is where the two
  packages compose: the Jacobian of a diagonal block is diagonal, which is
  exactly a scalar link's contract. `scaled_matrix()` is a fixed matrix
  carried by one scale, which covers the ridge, every single-smoothing-
  parameter spline penalty and the fully known case.

* Rank-deficient matrices are admitted, because a spline penalty is singular
  by construction and that is what makes it a penalty rather than a density.
  `param_logdet()` returns the log pseudo-determinant there, and its
  derivative under a scaled parameter is the rank -- the term without which the
  scale is not estimable.

* The rank and the null space are computed once at construction, from the
  components, by `param_null_basis()`. Reading a rank off an assembled sum is
  not scale invariant: on a tensor-product penalty of true rank 28 out of 32,
  counting eigenvalues above a relative tolerance gives 28 at equal scales, 26
  at a ratio of 1e8 and 16 at 1e12, while the null-space residual does not
  move. Membership is therefore tested against the null basis.

* `check_parameter()` runs nine checks, each against a route the
  implementation does not take, and reports a quantity that came from the base
  class as not checked rather than as passed.
