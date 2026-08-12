# parameters7 0.11.0

* The numerical fallbacks take their stencils from numericals7. The package
  imported it for the enumerations and used none of its differentiation:
  `mixed_stencil()` carried a transcribed table of the offsets and weights for
  orders one to four, which is a second copy of `fd_weights()`, and every
  fallback wrote its own central difference. The table is read from
  numericals7 now, `fd_step()` delegates rather than restating
  `eps^(1/(order+2))`, and `fd_along()` is the one place a difference along a
  free value is assembled. The rules are the same rules, so every number is
  unchanged; what moves is where they live.

# parameters7 0.10.0

* `sum_struct(components, link)`: a non-negative combination of fixed
  symmetric positive semidefinite matrices, which is the variance-components
  covariance and also the matrix a penalty with one smoothing parameter per
  component assembles. The value is linear in the weights, so a derivative
  component is zero unless every index names the same free value. The
  log-determinant is not separable and its derivatives come from the cyclic
  trace expansion
  \eqn{(-1)^{n-1}\sum_\sigma \mathrm{tr}(M^{-1}P_{\sigma(1)}\cdots
  M^{-1}P_{\sigma(n)})}, carried onto the free scale by a chain rule with a
  diagonal Jacobian.

  The orderings are counted WITH MULTIPLICITY. Deduplicating the ones that
  coincide because an index repeats makes a third derivative in one weight
  too small by a factor of two and a fourth by six, which is exactly what the
  first check of the expansion reported before the check itself was fixed.

  The rank is fixed at construction from the components stacked and
  individually normalized, the null space of a sum of positive semidefinite
  matrices being the intersection of theirs. Measured on second and first
  differences over five coefficients, the null residual stays at 1e-16 with
  the weights fourteen orders of magnitude apart, where a count of small
  eigenvalues of the assembled matrix would not.

  This completes the three compositions of `piano_parameters7.txt`, after
  `kron_identity()` in 0.8.0 and `block_diag()` and `dr_prod()` in 0.9.0.

# parameters7 0.9.0

* `block_diag(...)`: a block diagonal of distinct matrix parameters, each
  carrying its own stretch of the free vector. Where `kron_identity()`
  replicates one block, this composes several, which is what a model with
  more than one random-effect term needs. Every component whose indices do
  not all belong to one block is exactly zero, at every order and for the
  log-determinant as well as for the value, because the free values of one
  block do not enter another; the rank and the null basis come from the
  components rather than from an assembled matrix.

* `dr_prod(dimension, correlation, link)`: a covariance written as
  \eqn{D R D}, so that the coordinates are the standard deviations and the
  correlations rather than a function of them. Every derivative factorizes
  as \eqn{\partial^{S_D}(d_i d_j)\,\partial^{S_R} R_{ij}}, the two groups
  of free values being disjoint, and the log-determinant
  \eqn{2\sum_j \log d_j + \log|R|} is separable in the scales and separable
  from the correlation. The correlation block must have full rank: a
  deficient one would give the product the null space \eqn{D^{-1}\ker R},
  which moves with the free vector, while the class records the rank and the
  null space as properties of the family.

# parameters7 0.8.0

* kron_identity(structure, m): m identical diagonal blocks of one
  matrix parameter sharing a free vector, every contract quantity a
  linear lift of the inner one. The first composition wrapper; built
  for grouped random effects.

# parameters7 0.7.0

* The log-Cholesky derivative assembly is compiled. Every derivative of
  the factor is a single-entry matrix, so each Leibniz term is one row,
  one column or one cell rather than a dense product; `chol_leibniz_cpp`
  exploits that and serves all four orders, `param_d1` and `param_d2`
  included. Measured at p = 8: order 4 from 8.44 s to 0.28 s (30x, the
  remainder being the 82,251 result matrices the contract returns),
  order 2 from 6 ms to 1.2 ms. The dense R assembly stays as the twin
  `.chol_leibniz_r`, compared at machine precision in the tests.

# parameters7 0.6.0

* `autoregressive()` no longer consumes jets: the Levinson-Durbin recursion
  now propagates its derivative arrays in compiled code (`ar_taylor_cpp`),
  each tracked quantity holding its value and full symmetric tensors to
  fourth order, combined by the product rule written out per order. The
  package gains its first compiled code, and Rcpp enters Imports and
  LinkingTo. `ar_prediction()` and the log-determinant, which need no
  derivatives, read the link inverses directly. The q = 1 derivatives are
  pinned in the tests against products of link derivatives, formulas that
  share no code with the kernel.

# parameters7 0.5.0

* The jets move to numericals7, where the toolkit's numerical layer now
  lives: nothing about them concerned constrained parameters, and they landed
  here only because their first consumer -- the Levinson-Durbin recursion --
  did. `autoregressive()` now consumes them through numericals7, and its
  recursion is written on the arithmetic operators the jets dispatch on,
  which is what they exist for. A clean cut: `jet_layout()`, `jet_var()` and
  `set_partitions()` are no longer exported from here, and `tuple_indices()`
  delegates to the one enumeration in numericals7, so the two copies that
  could have disagreed no longer exist.

# parameters7 0.4.0

* Jets gain a composition rule and a vocabulary of functions, and can be
  written in ordinary R. `jet_compose()` applies a smooth scalar function to a
  jet given that function's own five derivatives, summing over the set
  partitions of the positions of an index tuple, so `exp`, `log`, an arbitrary
  power, `sqrt`, `gamma`, `lgamma`, `digamma` and `trigamma` each cost five
  numbers rather than a chain rule of their own. `Ops` and `Math` dispatch on
  the class, so a map reads `mu / gamma(1 + 1 / sigma)` and carries every
  partial derivative to fourth order with it. Comparison operators and the
  non-smooth functions are rejected: a branch taken on a jet would keep one
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
  correlation matrix in the spherical parametrization of Rapisarda, Brigo and
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
  in the parametrization of Pinheiro and Bates (1996), whose log-determinant
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
