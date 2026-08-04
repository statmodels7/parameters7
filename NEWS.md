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
