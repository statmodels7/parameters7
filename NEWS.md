# covstructs7 0.1.0

* First release. `covstruct` is the abstract class of constrained matrix
  parameters: a map from an unconstrained vector to a symmetric matrix,
  together with its first and second derivatives, its log-determinant and the
  solves a likelihood asks of it. Only `struct_matrix()` is compulsory;
  everything else has a method registered on the base class, so a new
  structure is a subclass and one method.

* Four families. `log_cholesky()` is the unstructured positive definite case
  in the parametrisation of Pinheiro and Bates (1996), whose log-determinant
  is linear in the free vector. `diag_struct()` and `scalar_struct()` carry
  their entries through a `linkfunctions7` link, which is where the two
  packages compose: the Jacobian of a diagonal block is diagonal, which is
  exactly a scalar link's contract. `scaled_struct()` is a fixed matrix
  carried by one scale, which covers the ridge, every single-smoothing-
  parameter spline penalty and the fully known case.

* Rank-deficient matrices are admitted, because a spline penalty is singular
  by construction and that is what makes it a penalty rather than a density.
  `struct_logdet()` returns the log pseudo-determinant there, and its
  derivative under a scaled structure is the rank -- the term without which the
  scale is not estimable.

* The rank and the null space are computed once at construction, from the
  components, by `struct_null_basis()`. Reading a rank off an assembled sum is
  not scale invariant: on a tensor-product penalty of true rank 28 out of 32,
  counting eigenvalues above a relative tolerance gives 28 at equal scales, 26
  at a ratio of 1e8 and 16 at 1e12, while the null-space residual does not
  move. Membership is therefore tested against the null basis.

* `check_covstruct()` runs nine checks, each against a route the
  implementation does not take, and reports a quantity that came from the base
  class as not checked rather than as passed.
