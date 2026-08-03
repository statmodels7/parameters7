# Package index

## Families

The parametrisations. Each maps an unconstrained vector to a matrix in
its own constrained set, and each owns the side of that matrix.

- [`log_cholesky()`](https://statmodels7.github.io/covstructs7/reference/log_cholesky.md)
  : Construct an Unstructured Positive Definite Structure
- [`diag_struct()`](https://statmodels7.github.io/covstructs7/reference/diag_struct.md)
  : Construct a Diagonal Structure
- [`scalar_struct()`](https://statmodels7.github.io/covstructs7/reference/scalar_struct.md)
  : Construct a Scalar Multiple of the Identity
- [`scaled_struct()`](https://statmodels7.github.io/covstructs7/reference/scaled_struct.md)
  : Construct a Scaled Fixed Matrix

## The map and its inverse

- [`struct_matrix()`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md)
  : The Matrix a Structure Produces
- [`struct_free()`](https://statmodels7.github.io/covstructs7/reference/struct_free.md)
  : The Free Vector Behind a Matrix

## Derivatives

Exact to second order, on the unconstrained scale, which is what joint
estimation over the coefficients and the matrix parameters consumes.

- [`struct_dmatrix()`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.md)
  : First Derivatives of a Structure's Matrix
- [`struct_d2matrix()`](https://statmodels7.github.io/covstructs7/reference/struct_d2matrix.md)
  : Second Derivatives of a Structure's Matrix
- [`struct_pair_names()`](https://statmodels7.github.io/covstructs7/reference/struct_pair_names.md)
  : Names of the Distinct Second-Derivative Components
- [`struct_pair_indices()`](https://statmodels7.github.io/covstructs7/reference/struct_pair_indices.md)
  : Index Pairs Behind the Second-Derivative Names

## What a likelihood asks of the matrix

The log-determinant, or the log pseudo-determinant when the family is
rank deficient, and the solves – computed through a factor rather than
through an explicit inverse.

- [`struct_logdet()`](https://statmodels7.github.io/covstructs7/reference/struct_logdet.md)
  : Log-Determinant of a Structure's Matrix
- [`struct_dlogdet()`](https://statmodels7.github.io/covstructs7/reference/struct_dlogdet.md)
  : Gradient of the Log-Determinant
- [`struct_d2logdet()`](https://statmodels7.github.io/covstructs7/reference/struct_d2logdet.md)
  : Hessian of the Log-Determinant
- [`struct_solve()`](https://statmodels7.github.io/covstructs7/reference/struct_solve.md)
  : Solve Through a Structure's Matrix
- [`struct_factor()`](https://statmodels7.github.io/covstructs7/reference/struct_factor.md)
  : A Factor of a Structure's Matrix

## Rank

A rank read off an assembled matrix is not scale invariant, so it is
computed once from the components and membership is tested against the
null space it implies.

- [`struct_null_basis()`](https://statmodels7.github.io/covstructs7/reference/struct_null_basis.md)
  : An Orthonormal Basis of a Null Space, and the Rank That Goes With It

## Tools

- [`check_covstruct()`](https://statmodels7.github.io/covstructs7/reference/check_covstruct.md)
  : Validate a Covariance Structure
- [`struct_is_numerical()`](https://statmodels7.github.io/covstructs7/reference/struct_is_numerical.md)
  : Which of a Structure's Quantities Come From the Base Class
- [`numerical_dmatrix()`](https://statmodels7.github.io/covstructs7/reference/numerical_dmatrix.md)
  : Numerical First Derivatives of a Structure's Matrix
- [`numerical_d2matrix()`](https://statmodels7.github.io/covstructs7/reference/numerical_d2matrix.md)
  : Numerical Second Derivatives of a Structure's Matrix

## Classes

The S7 classes; each page lists the methods that dispatch on it.

- [`covstruct()`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  : Constrained Matrix Parameter
- [`LogCholeskyStruct()`](https://statmodels7.github.io/covstructs7/reference/LogCholeskyStruct.md)
  : Unstructured Positive Definite Structure
- [`DiagStruct()`](https://statmodels7.github.io/covstructs7/reference/DiagStruct.md)
  : Diagonal Structure
- [`ScaledStruct()`](https://statmodels7.github.io/covstructs7/reference/ScaledStruct.md)
  : Scaled Fixed Matrix Structure

## Base-class defaults

The methods registered on the abstract class, which every structure
inherits unless it registers something more specific.

- [`check_covstruct()`](https://statmodels7.github.io/covstructs7/reference/check_covstruct.md)
  : Validate a Covariance Structure

## Internals

The machinery the exported functions are built from. None of it is
exported, and none is needed to use the package; it is documented so
that the derivations can be followed from the code that implements them.

- [`check_eta()`](https://statmodels7.github.io/covstructs7/reference/check_eta.md)
  : Validate a Free Vector Against a Structure
- [`check_matrix()`](https://statmodels7.github.io/covstructs7/reference/check_matrix.md)
  : Validate a Matrix Handed Back to a Structure
- [`check_positive_link()`](https://statmodels7.github.io/covstructs7/reference/check_positive_link.md)
  : Refuse a Link That Does Not Reach the Positive Half Line
- [`check_row()`](https://statmodels7.github.io/covstructs7/reference/check_row.md)
  : One Row of a Diagnostic Table
- [`check_struct_args()`](https://statmodels7.github.io/covstructs7/reference/check_struct_args.md)
  : Validate the Arguments Shared by Every Constructor
- [`chol_assemble()`](https://statmodels7.github.io/covstructs7/reference/chol_assemble.md)
  : The Cholesky Factor Behind a Free Vector
- [`chol_pd()`](https://statmodels7.github.io/covstructs7/reference/chol_pd.md)
  : Cholesky Factorisation, With the Rank Decided Before It
- [`chol_positions()`](https://statmodels7.github.io/covstructs7/reference/chol_positions.md)
  : Positions of the Free Values in the Lower Triangle
- [`covstructs7`](https://statmodels7.github.io/covstructs7/reference/covstructs7-package.md)
  [`covstructs7-package`](https://statmodels7.github.io/covstructs7/reference/covstructs7-package.md)
  : covstructs7: An S7 Framework for Constrained Matrix Parameters
- [`diag_entries()`](https://statmodels7.github.io/covstructs7/reference/diag_entries.md)
  : The Diagonal Entries Behind a Free Vector
- [`diag_multiplicity()`](https://statmodels7.github.io/covstructs7/reference/diag_multiplicity.md)
  : The Number of Diagonal Entries Each Free Value Owns
- [`diag_owner()`](https://statmodels7.github.io/covstructs7/reference/diag_owner.md)
  : The Free Value Each Diagonal Entry Belongs To
- [`empty_null_basis()`](https://statmodels7.github.io/covstructs7/reference/empty_null_basis.md)
  : No Null Space
- [`fd_step()`](https://statmodels7.github.io/covstructs7/reference/fd_step.md)
  : Finite-Difference Step for a Free Value
- [`is_base_struct_class()`](https://statmodels7.github.io/covstructs7/reference/is_base_struct_class.md)
  : Is This the Package's Own Base Class?
- [`name_dims()`](https://statmodels7.github.io/covstructs7/reference/name_dims.md)
  : Name the Rows and Columns of a Structure's Matrix
- [`print.covstruct`](https://statmodels7.github.io/covstructs7/reference/print.covstruct.md)
  : Print a Covariance Structure
- [`scaled_scale()`](https://statmodels7.github.io/covstructs7/reference/scaled_scale.md)
  : The Scale Behind a Free Vector, and Its Derivatives
- [`spectrum_pinv()`](https://statmodels7.github.io/covstructs7/reference/spectrum_pinv.md)
  : Moore-Penrose Inverse From a Structure's Spectrum
- [`struct_d2logdet.DiagStruct`](https://statmodels7.github.io/covstructs7/reference/struct_d2logdet.DiagStruct.md)
  : Log-Determinant Hessian of a Diagonal Structure
- [`struct_d2logdet.LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/struct_d2logdet.LogCholeskyStruct.md)
  : Log-Determinant Hessian of a Log-Cholesky Structure
- [`struct_d2logdet.ScaledStruct`](https://statmodels7.github.io/covstructs7/reference/struct_d2logdet.ScaledStruct.md)
  : Log-Determinant Hessian of a Scaled Structure
- [`struct_d2logdet.covstruct`](https://statmodels7.github.io/covstructs7/reference/struct_d2logdet.covstruct.md)
  : Default Log-Determinant Hessian
- [`struct_d2matrix.DiagStruct`](https://statmodels7.github.io/covstructs7/reference/struct_d2matrix.DiagStruct.md)
  : Second Derivatives of a Diagonal Structure
- [`struct_d2matrix.LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/struct_d2matrix.LogCholeskyStruct.md)
  : Second Derivatives of a Log-Cholesky Structure
- [`struct_d2matrix.ScaledStruct`](https://statmodels7.github.io/covstructs7/reference/struct_d2matrix.ScaledStruct.md)
  : Second Derivative of a Scaled Structure
- [`struct_d2matrix.covstruct`](https://statmodels7.github.io/covstructs7/reference/struct_d2matrix.covstruct.md)
  : Default Second Derivatives
- [`struct_dlogdet.DiagStruct`](https://statmodels7.github.io/covstructs7/reference/struct_dlogdet.DiagStruct.md)
  : Log-Determinant Gradient of a Diagonal Structure
- [`struct_dlogdet.LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/struct_dlogdet.LogCholeskyStruct.md)
  : Log-Determinant Gradient of a Log-Cholesky Structure
- [`struct_dlogdet.ScaledStruct`](https://statmodels7.github.io/covstructs7/reference/struct_dlogdet.ScaledStruct.md)
  : Log-Determinant Gradient of a Scaled Structure
- [`struct_dlogdet.covstruct`](https://statmodels7.github.io/covstructs7/reference/struct_dlogdet.covstruct.md)
  : Default Log-Determinant Gradient
- [`struct_dmatrix.DiagStruct`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.DiagStruct.md)
  : First Derivatives of a Diagonal Structure
- [`struct_dmatrix.LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.LogCholeskyStruct.md)
  : First Derivatives of a Log-Cholesky Structure
- [`struct_dmatrix.ScaledStruct`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.ScaledStruct.md)
  : First Derivative of a Scaled Structure
- [`struct_dmatrix.covstruct`](https://statmodels7.github.io/covstructs7/reference/struct_dmatrix.covstruct.md)
  : Default First Derivatives
- [`struct_factor.DiagStruct`](https://statmodels7.github.io/covstructs7/reference/struct_factor.DiagStruct.md)
  : Factor of a Diagonal Structure
- [`struct_factor.LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/struct_factor.LogCholeskyStruct.md)
  : Factor of a Log-Cholesky Structure
- [`struct_factor.covstruct`](https://statmodels7.github.io/covstructs7/reference/struct_factor.covstruct.md)
  : Default Factor
- [`struct_free.DiagStruct`](https://statmodels7.github.io/covstructs7/reference/struct_free.DiagStruct.md)
  : Free Vector of a Diagonal Structure
- [`struct_free.LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/struct_free.LogCholeskyStruct.md)
  : Free Vector of a Log-Cholesky Structure
- [`struct_free.ScaledStruct`](https://statmodels7.github.io/covstructs7/reference/struct_free.ScaledStruct.md)
  : Free Vector of a Scaled Structure
- [`struct_free.covstruct`](https://statmodels7.github.io/covstructs7/reference/struct_free.covstruct.md)
  : Refusal to Invert Without a Closed Form
- [`struct_logdet.DiagStruct`](https://statmodels7.github.io/covstructs7/reference/struct_logdet.DiagStruct.md)
  : Log-Determinant of a Diagonal Structure
- [`struct_logdet.LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/struct_logdet.LogCholeskyStruct.md)
  : Log-Determinant of a Log-Cholesky Structure
- [`struct_logdet.ScaledStruct`](https://statmodels7.github.io/covstructs7/reference/struct_logdet.ScaledStruct.md)
  : Log-Determinant of a Scaled Structure
- [`struct_logdet.covstruct`](https://statmodels7.github.io/covstructs7/reference/struct_logdet.covstruct.md)
  : Default Log-Determinant
- [`struct_matrix.DiagStruct`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.DiagStruct.md)
  : Matrix of a Diagonal Structure
- [`struct_matrix.LogCholeskyStruct`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.LogCholeskyStruct.md)
  : Matrix of a Log-Cholesky Structure
- [`struct_matrix.ScaledStruct`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.ScaledStruct.md)
  : Matrix of a Scaled Structure
- [`struct_solve.covstruct`](https://statmodels7.github.io/covstructs7/reference/struct_solve.covstruct.md)
  : Default Solve
- [`struct_spectrum()`](https://statmodels7.github.io/covstructs7/reference/struct_spectrum.md)
  : The Spectral Decomposition a Structure's Quantities Are Read From
- [`sweep_etas()`](https://statmodels7.github.io/covstructs7/reference/sweep_etas.md)
  : Free Vectors to Sweep a Structure Over
