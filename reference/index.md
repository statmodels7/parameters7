# Package index

## Families

The parametrisations. Each maps an unconstrained vector to a matrix in
its own constrained set, and each owns the side of that matrix.

- [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  : Construct an Unstructured Positive Definite Parameter
- [`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
  : Construct a Matrix Logarithm Parameter
- [`diagonal_matrix()`](https://statmodels7.github.io/parameters7/reference/diagonal_matrix.md)
  : Construct a Diagonal Parameter
- [`scalar_matrix()`](https://statmodels7.github.io/parameters7/reference/scalar_matrix.md)
  : Construct a Scalar Multiple of the Identity
- [`scaled_matrix()`](https://statmodels7.github.io/parameters7/reference/scaled_matrix.md)
  : Construct a Scaled Fixed Matrix
- [`simplex()`](https://statmodels7.github.io/parameters7/reference/simplex.md)
  : Construct a Simplex Parameter
- [`transition_matrix()`](https://statmodels7.github.io/parameters7/reference/transition_matrix.md)
  : Construct a Transition Matrix Parameter

## The map and its inverse

- [`param_value()`](https://statmodels7.github.io/parameters7/reference/param_value.md)
  : The Matrix a Parameter Produces
- [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
  : The Free Vector Behind a Matrix

## Derivatives

Exact to fourth order, on the unconstrained scale, which is what joint
estimation over the coefficients and the constrained parameters
consumes.

- [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
  : First Derivatives of a Parameter's Matrix
- [`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
  : Second Derivatives of a Parameter's Matrix
- [`param_d3()`](https://statmodels7.github.io/parameters7/reference/param_d3.md)
  : Third Derivatives of a Parameter's Value
- [`param_d4()`](https://statmodels7.github.io/parameters7/reference/param_d4.md)
  : Fourth Derivatives of a Parameter's Value
- [`param_tuple_names()`](https://statmodels7.github.io/parameters7/reference/param_tuple_names.md)
  : Names of the Distinct Derivative Components
- [`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
  : Index Tuples Behind the Derivative Component Names

## What a likelihood asks of the matrix

The log-determinant, or the log pseudo-determinant when the family is
rank deficient, and the solves – computed through a factor rather than
through an explicit inverse.

- [`param_logdet()`](https://statmodels7.github.io/parameters7/reference/param_logdet.md)
  : Log-Determinant of a Parameter's Matrix
- [`param_dlogdet()`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.md)
  : Gradient of the Log-Determinant
- [`param_d2logdet()`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.md)
  : Hessian of the Log-Determinant
- [`param_d3logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
  [`param_d4logdet()`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.md)
  : Third and Fourth Derivatives of the Log-Determinant
- [`param_solve()`](https://statmodels7.github.io/parameters7/reference/param_solve.md)
  : Solve Through a Parameter's Matrix
- [`param_factor()`](https://statmodels7.github.io/parameters7/reference/param_factor.md)
  : A Factor of a Parameter's Matrix

## Rank

A rank read off an assembled matrix is not scale invariant, so it is
computed once from the components and membership is tested against the
null space it implies.

- [`param_null_basis()`](https://statmodels7.github.io/parameters7/reference/param_null_basis.md)
  : An Orthonormal Basis of a Null Space, and the Rank That Goes With It

## Tools

- [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  : Validate a Covariance Parameter
- [`param_is_numerical()`](https://statmodels7.github.io/parameters7/reference/param_is_numerical.md)
  : Which of a Parameter's Quantities Come From the Base Class
- [`numerical_d1()`](https://statmodels7.github.io/parameters7/reference/numerical_d1.md)
  : Numerical First Derivatives of a Parameter's Matrix
- [`numerical_d2()`](https://statmodels7.github.io/parameters7/reference/numerical_d2.md)
  : Numerical Second Derivatives of a Parameter's Matrix
- [`numerical_d3()`](https://statmodels7.github.io/parameters7/reference/numerical_d3.md)
  : Numerical Third Derivatives of a Parameter's Value
- [`numerical_d4()`](https://statmodels7.github.io/parameters7/reference/numerical_d4.md)
  : Numerical Fourth Derivatives of a Parameter's Value

## Classes

The S7 classes; each page lists the methods that dispatch on it.

- [`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
  : Constrained Matrix Parameter
- [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  : Constrained Symmetric Matrix Parameter
- [`LogCholeskyParam()`](https://statmodels7.github.io/parameters7/reference/LogCholeskyParam.md)
  : Unstructured Positive Definite Parameter
- [`MatrixLogParam()`](https://statmodels7.github.io/parameters7/reference/MatrixLogParam.md)
  : Matrix Logarithm Parameter
- [`SimplexParam()`](https://statmodels7.github.io/parameters7/reference/SimplexParam.md)
  : Simplex Parameter
- [`TransitionMatrixParam()`](https://statmodels7.github.io/parameters7/reference/TransitionMatrixParam.md)
  : Transition Matrix Parameter
- [`DiagMatrixParam()`](https://statmodels7.github.io/parameters7/reference/DiagMatrixParam.md)
  : Diagonal Parameter
- [`ScaledMatrixParam()`](https://statmodels7.github.io/parameters7/reference/ScaledMatrixParam.md)
  : Scaled Fixed Matrix Parameter

## Base-class defaults

The methods registered on the abstract class, which every parameter
inherits unless it registers something more specific.

- [`check_parameter()`](https://statmodels7.github.io/parameters7/reference/check_parameter.md)
  : Validate a Covariance Parameter
- [`matrix_parameter()`](https://statmodels7.github.io/parameters7/reference/matrix_parameter.md)
  : Constrained Symmetric Matrix Parameter

## Internals

The machinery the exported functions are built from. None of it is
exported, and none is needed to use the package; it is documented so
that the derivations can be followed from the code that implements them.

- [`check_eta()`](https://statmodels7.github.io/parameters7/reference/check_eta.md)
  : Validate a Free Vector Against a Parameter
- [`check_matrix()`](https://statmodels7.github.io/parameters7/reference/check_matrix.md)
  : Validate a Matrix Handed Back to a Parameter
- [`check_param_args()`](https://statmodels7.github.io/parameters7/reference/check_param_args.md)
  : Validate the Arguments Shared by Every Constructor
- [`check_parameter_vector()`](https://statmodels7.github.io/parameters7/reference/check_parameter_vector.md)
  : The Reduced Battery for a Parameter That Is Not a Matrix
- [`check_positive_link()`](https://statmodels7.github.io/parameters7/reference/check_positive_link.md)
  : Refuse a Link That Does Not Reach the Positive Half Line
- [`check_row()`](https://statmodels7.github.io/parameters7/reference/check_row.md)
  : One Row of a Diagnostic Table
- [`chol_assemble()`](https://statmodels7.github.io/parameters7/reference/chol_assemble.md)
  : The Cholesky Factor Behind a Free Vector
- [`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)
  : Cholesky Factorisation, With the Rank Decided Before It
- [`chol_positions()`](https://statmodels7.github.io/parameters7/reference/chol_positions.md)
  : Positions of the Free Values in the Lower Triangle
- [`combinat_perms()`](https://statmodels7.github.io/parameters7/reference/combinat_perms.md)
  : All Orderings of a Tuple
- [`dd_exp()`](https://statmodels7.github.io/parameters7/reference/dd_exp.md)
  : Divided Differences of the Exponential
- [`diag_entries()`](https://statmodels7.github.io/parameters7/reference/diag_entries.md)
  : The Diagonal Entries Behind a Free Vector
- [`diag_multiplicity()`](https://statmodels7.github.io/parameters7/reference/diag_multiplicity.md)
  : The Number of Diagonal Entries Each Free Value Owns
- [`diag_owner()`](https://statmodels7.github.io/parameters7/reference/diag_owner.md)
  : The Free Value Each Diagonal Entry Belongs To
- [`empty_null_basis()`](https://statmodels7.github.io/parameters7/reference/empty_null_basis.md)
  : No Null Space
- [`fd_step()`](https://statmodels7.github.io/parameters7/reference/fd_step.md)
  : Finite-Difference Step for a Free Value
- [`is_base_param_class()`](https://statmodels7.github.io/parameters7/reference/is_base_param_class.md)
  : Is This the Package's Own Base Class?
- [`mixed_stencil()`](https://statmodels7.github.io/parameters7/reference/mixed_stencil.md)
  : One Product Stencil for a Mixed Partial Derivative
- [`mlog_basis()`](https://statmodels7.github.io/parameters7/reference/mlog_basis.md)
  : The Basis Direction of One Free Value
- [`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md)
  : One Derivative Component by the Daleckii-Krein Contraction
- [`mlog_expm_small()`](https://statmodels7.github.io/parameters7/reference/mlog_expm_small.md)
  : Exponential of a Small Upper Triangular Matrix
- [`mlog_s()`](https://statmodels7.github.io/parameters7/reference/mlog_s.md)
  : The Symmetric Matrix Behind a Free Vector
- [`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md)
  : The Eigendecomposition and Divided-Difference Tables at a Point
- [`name_dims()`](https://statmodels7.github.io/parameters7/reference/name_dims.md)
  : Name the Rows and Columns of a Parameter's Matrix
- [`param_d1.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d1.DiagMatrixParam.md)
  : First Derivatives of a Diagonal Parameter
- [`param_d1.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_d1.LogCholeskyParam.md)
  : First Derivatives of a Log-Cholesky Parameter
- [`param_d1.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_d1.MatrixLogParam.md)
  : First Derivatives of a Matrix Logarithm Parameter
- [`param_d1.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d1.ScaledMatrixParam.md)
  : First Derivative of a Scaled Parameter
- [`param_d1.SimplexParam`](https://statmodels7.github.io/parameters7/reference/param_d1.SimplexParam.md)
  : First Derivatives of a Simplex Parameter
- [`param_d1.TransitionMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d1.TransitionMatrixParam.md)
  : First Derivatives of a Transition Matrix Parameter
- [`param_d1.parameter`](https://statmodels7.github.io/parameters7/reference/param_d1.parameter.md)
  : Default First Derivatives
- [`param_d2.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d2.DiagMatrixParam.md)
  : Second Derivatives of a Diagonal Parameter
- [`param_d2.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_d2.LogCholeskyParam.md)
  : Second Derivatives of a Log-Cholesky Parameter
- [`param_d2.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_d2.MatrixLogParam.md)
  : Second Derivatives of a Matrix Logarithm Parameter
- [`param_d2.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d2.ScaledMatrixParam.md)
  : Second Derivative of a Scaled Parameter
- [`param_d2.SimplexParam`](https://statmodels7.github.io/parameters7/reference/param_d2.SimplexParam.md)
  : Second Derivatives of a Simplex Parameter
- [`param_d2.TransitionMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d2.TransitionMatrixParam.md)
  : Second Derivatives of a Transition Matrix Parameter
- [`param_d2.parameter`](https://statmodels7.github.io/parameters7/reference/param_d2.parameter.md)
  : Default Second Derivatives
- [`param_d2logdet.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.DiagMatrixParam.md)
  : Log-Determinant Hessian of a Diagonal Parameter
- [`param_d2logdet.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.LogCholeskyParam.md)
  : Log-Determinant Hessian of a Log-Cholesky Parameter
- [`param_d2logdet.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.MatrixLogParam.md)
  : Log-Determinant Hessian of a Matrix Logarithm Parameter
- [`param_d2logdet.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.ScaledMatrixParam.md)
  : Log-Determinant Hessian of a Scaled Parameter
- [`param_d2logdet.parameter`](https://statmodels7.github.io/parameters7/reference/param_d2logdet.parameter.md)
  : Default Log-Determinant Hessian
- [`param_d3.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d3.DiagMatrixParam.md)
  : Third Derivatives of a Diagonal Parameter
- [`param_d3.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_d3.LogCholeskyParam.md)
  : Third Derivatives of a Log-Cholesky Parameter
- [`param_d3.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_d3.MatrixLogParam.md)
  : Third Derivatives of a Matrix Logarithm Parameter
- [`param_d3.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d3.ScaledMatrixParam.md)
  : Third Derivatives of a Scaled Parameter
- [`param_d3.SimplexParam`](https://statmodels7.github.io/parameters7/reference/param_d3.SimplexParam.md)
  : Third Derivatives of a Simplex Parameter
- [`param_d3.TransitionMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d3.TransitionMatrixParam.md)
  : Third Derivatives of a Transition Matrix Parameter
- [`param_d3.parameter`](https://statmodels7.github.io/parameters7/reference/param_d3.parameter.md)
  : Default Third Derivatives
- [`param_d3logdet.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.DiagMatrixParam.md)
  : Higher Log-Determinant Derivatives of a Diagonal Parameter
- [`param_d3logdet.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.LogCholeskyParam.md)
  : Higher Log-Determinant Derivatives of a Log-Cholesky Parameter
- [`param_d3logdet.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.MatrixLogParam.md)
  : Higher Log-Determinant Derivatives of a Matrix Logarithm Parameter
- [`param_d3logdet.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.ScaledMatrixParam.md)
  : Higher Log-Determinant Derivatives of a Scaled Parameter
- [`param_d3logdet.matrix_parameter`](https://statmodels7.github.io/parameters7/reference/param_d3logdet.matrix_parameter.md)
  : Default Higher Log-Determinant Derivatives
- [`param_d4.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d4.DiagMatrixParam.md)
  : Fourth Derivatives of a Diagonal Parameter
- [`param_d4.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_d4.LogCholeskyParam.md)
  : Fourth Derivatives of a Log-Cholesky Parameter
- [`param_d4.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_d4.MatrixLogParam.md)
  : Fourth Derivatives of a Matrix Logarithm Parameter
- [`param_d4.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d4.ScaledMatrixParam.md)
  : Fourth Derivatives of a Scaled Parameter
- [`param_d4.SimplexParam`](https://statmodels7.github.io/parameters7/reference/param_d4.SimplexParam.md)
  : Fourth Derivatives of a Simplex Parameter
- [`param_d4.TransitionMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d4.TransitionMatrixParam.md)
  : Fourth Derivatives of a Transition Matrix Parameter
- [`param_d4.parameter`](https://statmodels7.github.io/parameters7/reference/param_d4.parameter.md)
  : Default Fourth Derivatives
- [`param_d4logdet.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.DiagMatrixParam.md)
  : Fourth Log-Determinant Derivatives of a Diagonal Parameter
- [`param_d4logdet.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.LogCholeskyParam.md)
  : Fourth Log-Determinant Derivatives of a Log-Cholesky Parameter
- [`param_d4logdet.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.MatrixLogParam.md)
  : Fourth Log-Determinant Derivatives of a Matrix Logarithm Parameter
- [`param_d4logdet.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.ScaledMatrixParam.md)
  : Fourth Log-Determinant Derivatives of a Scaled Parameter
- [`param_d4logdet.matrix_parameter`](https://statmodels7.github.io/parameters7/reference/param_d4logdet.matrix_parameter.md)
  : Default Fourth Log-Determinant Derivatives
- [`param_dlogdet.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.DiagMatrixParam.md)
  : Log-Determinant Gradient of a Diagonal Parameter
- [`param_dlogdet.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.LogCholeskyParam.md)
  : Log-Determinant Gradient of a Log-Cholesky Parameter
- [`param_dlogdet.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.MatrixLogParam.md)
  : Log-Determinant Gradient of a Matrix Logarithm Parameter
- [`param_dlogdet.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.ScaledMatrixParam.md)
  : Log-Determinant Gradient of a Scaled Parameter
- [`param_dlogdet.parameter`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.parameter.md)
  : Default Log-Determinant Gradient
- [`param_factor.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_factor.DiagMatrixParam.md)
  : Factor of a Diagonal Parameter
- [`param_factor.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_factor.LogCholeskyParam.md)
  : Factor of a Log-Cholesky Parameter
- [`param_factor.parameter`](https://statmodels7.github.io/parameters7/reference/param_factor.parameter.md)
  : Default Factor
- [`param_free.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_free.DiagMatrixParam.md)
  : Free Vector of a Diagonal Parameter
- [`param_free.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_free.LogCholeskyParam.md)
  : Free Vector of a Log-Cholesky Parameter
- [`param_free.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_free.MatrixLogParam.md)
  : Free Vector of a Matrix Logarithm Parameter
- [`param_free.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_free.ScaledMatrixParam.md)
  : Free Vector of a Scaled Parameter
- [`param_free.SimplexParam`](https://statmodels7.github.io/parameters7/reference/param_free.SimplexParam.md)
  : Free Vector of a Simplex Parameter
- [`param_free.TransitionMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_free.TransitionMatrixParam.md)
  : Free Vector of a Transition Matrix Parameter
- [`param_free.parameter`](https://statmodels7.github.io/parameters7/reference/param_free.parameter.md)
  : Refusal to Invert Without a Closed Form
- [`param_logdet.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_logdet.DiagMatrixParam.md)
  : Log-Determinant of a Diagonal Parameter
- [`param_logdet.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_logdet.LogCholeskyParam.md)
  : Log-Determinant of a Log-Cholesky Parameter
- [`param_logdet.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_logdet.MatrixLogParam.md)
  : Log-Determinant of a Matrix Logarithm Parameter
- [`param_logdet.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_logdet.ScaledMatrixParam.md)
  : Log-Determinant of a Scaled Parameter
- [`param_logdet.parameter`](https://statmodels7.github.io/parameters7/reference/param_logdet.parameter.md)
  : Default Log-Determinant
- [`param_solve.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_solve.MatrixLogParam.md)
  : Solve of a Matrix Logarithm Parameter
- [`param_solve.parameter`](https://statmodels7.github.io/parameters7/reference/param_solve.parameter.md)
  : Default Solve
- [`param_spectrum()`](https://statmodels7.github.io/parameters7/reference/param_spectrum.md)
  : The Spectral Decomposition a Parameter's Quantities Are Read From
- [`param_value.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_value.DiagMatrixParam.md)
  : Matrix of a Diagonal Parameter
- [`param_value.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_value.LogCholeskyParam.md)
  : Matrix of a Log-Cholesky Parameter
- [`param_value.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_value.MatrixLogParam.md)
  : Value of a Matrix Logarithm Parameter
- [`param_value.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_value.ScaledMatrixParam.md)
  : Matrix of a Scaled Parameter
- [`param_value.SimplexParam`](https://statmodels7.github.io/parameters7/reference/param_value.SimplexParam.md)
  : Value of a Simplex Parameter
- [`param_value.TransitionMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_value.TransitionMatrixParam.md)
  : Value of a Transition Matrix Parameter
- [`parameters7`](https://statmodels7.github.io/parameters7/reference/parameters7-package.md)
  [`parameters7-package`](https://statmodels7.github.io/parameters7/reference/parameters7-package.md)
  : parameters7: An S7 Framework for Constrained Parameters
- [`print.parameter`](https://statmodels7.github.io/parameters7/reference/print.parameter.md)
  : Print a Covariance Parameter
- [`scaled_scale()`](https://statmodels7.github.io/parameters7/reference/scaled_scale.md)
  : The Scale Behind a Free Vector, and Its Derivatives
- [`simplex_components()`](https://statmodels7.github.io/parameters7/reference/simplex_components.md)
  : Extract Named Components From Softmax Tensors
- [`simplex_point()`](https://statmodels7.github.io/parameters7/reference/simplex_point.md)
  : The Softmax Point Behind a Free Vector
- [`simplex_tensors()`](https://statmodels7.github.io/parameters7/reference/simplex_tensors.md)
  : Derivative Tensors of the Softmax Map
- [`spectrum_pinv()`](https://statmodels7.github.io/parameters7/reference/spectrum_pinv.md)
  : Moore-Penrose Inverse From a Parameter's Spectrum
- [`sweep_etas()`](https://statmodels7.github.io/parameters7/reference/sweep_etas.md)
  : Free Vectors to Sweep a Parameter Over
- [`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md)
  : Derivative Components of a Transition Matrix Parameter
- [`tm_positions()`](https://statmodels7.github.io/parameters7/reference/tm_positions.md)
  : Row and Chart Coordinate of Each Free Value
