# Package index

## Families

The parametrizations. Each maps an unconstrained vector to a matrix in
its own constrained set, and each owns the side of that matrix.

- [`log_cholesky()`](https://statmodels7.github.io/parameters7/reference/log_cholesky.md)
  : Construct an Unstructured Positive Definite Parameter
- [`matrix_log()`](https://statmodels7.github.io/parameters7/reference/matrix_log.md)
  : Construct a Matrix Logarithm Parameter
- [`correlation_matrix()`](https://statmodels7.github.io/parameters7/reference/correlation_matrix.md)
  : Construct a Correlation Matrix Parameter
- [`compound_symmetry()`](https://statmodels7.github.io/parameters7/reference/compound_symmetry.md)
  : Construct a Compound Symmetry Parameter
- [`ar1()`](https://statmodels7.github.io/parameters7/reference/ar1.md)
  : Construct an AR(1) Parameter
- [`autoregressive()`](https://statmodels7.github.io/parameters7/reference/autoregressive.md)
  : Construct an Autoregressive Parameter
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
  : The Value a Parameter Produces
- [`param_free()`](https://statmodels7.github.io/parameters7/reference/param_free.md)
  : The Free Vector Behind a Value

## Derivatives

Exact to fourth order, on the unconstrained scale, which is what joint
estimation over the coefficients and the constrained parameters
consumes.

- [`param_d1()`](https://statmodels7.github.io/parameters7/reference/param_d1.md)
  : First Derivatives of a Parameter's Value
- [`param_d2()`](https://statmodels7.github.io/parameters7/reference/param_d2.md)
  : Second Derivatives of a Parameter's Value
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

## What a reader reads

The quantities a family is about, with the Jacobian of the map from the
free vector and the scale each interval belongs on, so that a consumer
can report them by the delta method rather than reporting coordinates.

- [`param_readable()`](https://statmodels7.github.io/parameters7/reference/param_readable.md)
  : Quantities a Family Is About

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
- [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  : Correlation Matrix Parameter
- [`CompoundSymmetryParam()`](https://statmodels7.github.io/parameters7/reference/CompoundSymmetryParam.md)
  : Compound Symmetry Parameter
- [`Ar1Param()`](https://statmodels7.github.io/parameters7/reference/Ar1Param.md)
  : AR(1) Parameter
- [`AutoregressiveParam()`](https://statmodels7.github.io/parameters7/reference/AutoregressiveParam.md)
  : Autoregressive Parameter
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

- [`ar1_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/ar1_logdet_terms.md)
  : The Log-Determinant Terms of an AR(1) Parameter
- [`ar1_pattern()`](https://statmodels7.github.io/parameters7/reference/ar1_pattern.md)
  : The Pattern of an AR(1) Parameter
- [`ar_assemble()`](https://statmodels7.github.io/parameters7/reference/ar_assemble.md)
  : The Matrix and Its Derivatives, From the Packed Arrays
- [`ar_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_derivative.md)
  : Derivative Components of an Autoregressive Parameter
- [`ar_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_logdet_derivative.md)
  : Log-Determinant Components of an Autoregressive Parameter
- [`ar_pack_col()`](https://statmodels7.github.io/parameters7/reference/ar_pack_col.md)
  : The Column of a Packed Derivative Record
- [`ar_prediction()`](https://statmodels7.github.io/parameters7/reference/ar_prediction.md)
  : The Prediction Form of an Autoregressive Parameter
- [`ar_taylor()`](https://statmodels7.github.io/parameters7/reference/ar_taylor.md)
  : The Levinson-Durbin Recursion With Its Derivatives
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
- [`chol_dfactor()`](https://statmodels7.github.io/parameters7/reference/chol_dfactor.md)
  : The Factor of a Log-Cholesky Parameter, and Its Derivatives
- [`chol_leibniz()`](https://statmodels7.github.io/parameters7/reference/chol_leibniz.md)
  : Derivative Components of a Log-Cholesky Parameter
- [`chol_pd()`](https://statmodels7.github.io/parameters7/reference/chol_pd.md)
  : Cholesky Factorization, With the Rank Decided Before It
- [`chol_positions()`](https://statmodels7.github.io/parameters7/reference/chol_positions.md)
  : Positions of the Free Values in the Lower Triangle
- [`combinat_perms()`](https://statmodels7.github.io/parameters7/reference/combinat_perms.md)
  : All Orderings of a Tuple
- [`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
  : Compose Two Scalar Maps, to Fourth Order
- [`corr_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_derivative.md)
  : Derivative Components of a Correlation Parameter
- [`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)
  : The Factor of a Correlation Parameter, and Its Derivatives
- [`corr_logdet_chains()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_chains.md)
  : Log-Determinant Chains of a Correlation Parameter
- [`corr_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_derivative.md)
  : Log-Determinant Components of a Correlation Parameter
- [`corr_tables()`](https://statmodels7.github.io/parameters7/reference/corr_tables.md)
  : Sines and Cosines of a Correlation Parameter's Angles
- [`cs_logdet_terms()`](https://statmodels7.github.io/parameters7/reference/cs_logdet_terms.md)
  : The Log-Determinant Terms of a Compound-Symmetric Parameter
- [`cs_pattern()`](https://statmodels7.github.io/parameters7/reference/cs_pattern.md)
  : The Pattern of a Compound-Symmetric Parameter
- [`dd_exp()`](https://statmodels7.github.io/parameters7/reference/dd_exp.md)
  : Divided Differences of the Exponential
- [`diag_dlog()`](https://statmodels7.github.io/parameters7/reference/diag_dlog.md)
  : Derivatives of the Logarithm of a Link
- [`diag_entries()`](https://statmodels7.github.io/parameters7/reference/diag_entries.md)
  : The Diagonal Entries Behind a Free Vector
- [`diag_higher()`](https://statmodels7.github.io/parameters7/reference/diag_higher.md)
  : Third and Fourth Derivatives of a Diagonal Matrix
- [`diag_logdet_higher()`](https://statmodels7.github.io/parameters7/reference/diag_logdet_higher.md)
  : Higher Derivatives of a Diagonal Log-Determinant
- [`diag_multiplicity()`](https://statmodels7.github.io/parameters7/reference/diag_multiplicity.md)
  : The Number of Diagonal Entries Each Free Value Owns
- [`diag_owner()`](https://statmodels7.github.io/parameters7/reference/diag_owner.md)
  : The Free Value Each Diagonal Entry Belongs To
- [`econ_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_derivative.md)
  : Derivative Components of an Economical Parameter
- [`econ_logdet_derivative()`](https://statmodels7.github.io/parameters7/reference/econ_logdet_derivative.md)
  : Log-Determinant Components of an Economical Parameter
- [`econ_scalars()`](https://statmodels7.github.io/parameters7/reference/econ_scalars.md)
  : The Scale and the Correlation of an Economical Parameter
- [`empty_null_basis()`](https://statmodels7.github.io/parameters7/reference/empty_null_basis.md)
  : No Null Space
- [`fd_step()`](https://statmodels7.github.io/parameters7/reference/fd_step.md)
  : Finite-Difference Step for a Free Value
- [`is_base_param_class()`](https://statmodels7.github.io/parameters7/reference/is_base_param_class.md)
  : Is This the Package's Own Base Class?
- [`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md)
  : A Gram Product's Derivatives From Its Factor's
- [`link_tag()`](https://statmodels7.github.io/parameters7/reference/link_tag.md)
  : A Short Name for a Link
- [`log_affine_derivs()`](https://statmodels7.github.io/parameters7/reference/log_affine_derivs.md)
  : Derivatives of a Sum of Logarithms of Affine Functions
- [`mixed_stencil()`](https://statmodels7.github.io/parameters7/reference/mixed_stencil.md)
  : One Product Stencil for a Mixed Partial Derivative
- [`mlog_basis()`](https://statmodels7.github.io/parameters7/reference/mlog_basis.md)
  : The Basis Direction of One Free Value
- [`mlog_contract()`](https://statmodels7.github.io/parameters7/reference/mlog_contract.md)
  : One Derivative Component by the Daleckii-Krein Contraction
- [`mlog_expm_small()`](https://statmodels7.github.io/parameters7/reference/mlog_expm_small.md)
  : Exponential of a Small Upper Triangular Matrix
- [`mlog_higher()`](https://statmodels7.github.io/parameters7/reference/mlog_higher.md)
  : Third and Fourth Derivatives of a Matrix Exponential
- [`mlog_s()`](https://statmodels7.github.io/parameters7/reference/mlog_s.md)
  : The Symmetric Matrix Behind a Free Vector
- [`mlog_tables()`](https://statmodels7.github.io/parameters7/reference/mlog_tables.md)
  : The Eigendecomposition and Divided-Difference Tables at a Point
- [`name_dims()`](https://statmodels7.github.io/parameters7/reference/name_dims.md)
  : Name the Rows and Columns of a Parameter's Matrix
- [`param_d1.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_d1.Ar1Param.md)
  [`param_d2.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_d1.Ar1Param.md)
  [`param_d3.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_d1.Ar1Param.md)
  [`param_d4.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_d1.Ar1Param.md)
  : Derivatives of an AR(1) Parameter
- [`param_d1.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_d1.AutoregressiveParam.md)
  [`param_d2.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_d1.AutoregressiveParam.md)
  [`param_d3.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_d1.AutoregressiveParam.md)
  [`param_d4.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_d1.AutoregressiveParam.md)
  : Derivatives of an Autoregressive Parameter
- [`param_d1.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_d1.CompoundSymmetryParam.md)
  [`param_d2.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_d1.CompoundSymmetryParam.md)
  [`param_d3.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_d1.CompoundSymmetryParam.md)
  [`param_d4.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_d1.CompoundSymmetryParam.md)
  : Derivatives of a Compound Symmetry Parameter
- [`param_d1.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_d1.CorrelationParam.md)
  : First Derivatives of a Correlation Parameter
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
- [`param_d2.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_d2.CorrelationParam.md)
  : Second Derivatives of a Correlation Parameter
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
- [`param_d3.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_d3.CorrelationParam.md)
  : Third Derivatives of a Correlation Parameter
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
- [`param_d4.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_d4.CorrelationParam.md)
  : Fourth Derivatives of a Correlation Parameter
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
- [`param_dlogdet.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.Ar1Param.md)
  [`param_d2logdet.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.Ar1Param.md)
  [`param_d3logdet.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.Ar1Param.md)
  [`param_d4logdet.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.Ar1Param.md)
  : Log-Determinant Derivatives of an AR(1) Parameter
- [`param_dlogdet.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md)
  [`param_d2logdet.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md)
  [`param_d3logdet.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md)
  [`param_d4logdet.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.AutoregressiveParam.md)
  : Log-Determinant Derivatives of an Autoregressive Parameter
- [`param_dlogdet.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CompoundSymmetryParam.md)
  [`param_d2logdet.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CompoundSymmetryParam.md)
  [`param_d3logdet.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CompoundSymmetryParam.md)
  [`param_d4logdet.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CompoundSymmetryParam.md)
  : Log-Determinant Derivatives of a Compound Symmetry Parameter
- [`param_dlogdet.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CorrelationParam.md)
  [`param_d2logdet.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CorrelationParam.md)
  [`param_d3logdet.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CorrelationParam.md)
  [`param_d4logdet.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_dlogdet.CorrelationParam.md)
  : Log-Determinant Derivatives of a Correlation Parameter
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
- [`param_factor.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_factor.CorrelationParam.md)
  : Factor of a Correlation Parameter
- [`param_factor.DiagMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_factor.DiagMatrixParam.md)
  : Factor of a Diagonal Parameter
- [`param_factor.LogCholeskyParam`](https://statmodels7.github.io/parameters7/reference/param_factor.LogCholeskyParam.md)
  : Factor of a Log-Cholesky Parameter
- [`param_factor.parameter`](https://statmodels7.github.io/parameters7/reference/param_factor.parameter.md)
  : Default Factor
- [`param_free.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_free.Ar1Param.md)
  : Free Vector of an AR(1) Parameter
- [`param_free.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_free.AutoregressiveParam.md)
  : Free Vector of an Autoregressive Parameter
- [`param_free.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_free.CompoundSymmetryParam.md)
  : Free Vector of a Compound Symmetry Parameter
- [`param_free.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_free.CorrelationParam.md)
  : Free Vector of a Correlation Parameter
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
- [`param_logdet.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_logdet.Ar1Param.md)
  : Log-Determinant of an AR(1) Parameter
- [`param_logdet.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_logdet.AutoregressiveParam.md)
  : Log-Determinant of an Autoregressive Parameter
- [`param_logdet.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_logdet.CompoundSymmetryParam.md)
  : Log-Determinant of a Compound Symmetry Parameter
- [`param_logdet.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_logdet.CorrelationParam.md)
  : Log-Determinant of a Correlation Parameter
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
- [`param_readable.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_readable.Ar1Param.md)
  : The Scale and the Correlation of an AR(1)
- [`param_readable.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_readable.AutoregressiveParam.md)
  : The Scale, the Partial Autocorrelations and the Coefficients
- [`param_readable.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_readable.CompoundSymmetryParam.md)
  : The Scale and the Common Correlation of a Compound Symmetry
- [`param_readable.ScaledMatrixParam`](https://statmodels7.github.io/parameters7/reference/param_readable.ScaledMatrixParam.md)
  : The Scale of a Fixed Matrix
- [`param_readable.SimplexParam`](https://statmodels7.github.io/parameters7/reference/param_readable.SimplexParam.md)
  : The Probabilities Behind an Additive Log-Ratio
- [`param_readable.parameter`](https://statmodels7.github.io/parameters7/reference/param_readable.parameter.md)
  : No Declared Quantities
- [`param_solve.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_solve.Ar1Param.md)
  : Solve of an AR(1) Parameter
- [`param_solve.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_solve.AutoregressiveParam.md)
  : Solve of an Autoregressive Parameter
- [`param_solve.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_solve.CompoundSymmetryParam.md)
  : Solve of a Compound Symmetry Parameter
- [`param_solve.MatrixLogParam`](https://statmodels7.github.io/parameters7/reference/param_solve.MatrixLogParam.md)
  : Solve of a Matrix Logarithm Parameter
- [`param_solve.parameter`](https://statmodels7.github.io/parameters7/reference/param_solve.parameter.md)
  : Default Solve
- [`param_spectrum()`](https://statmodels7.github.io/parameters7/reference/param_spectrum.md)
  : The Spectral Decomposition a Parameter's Quantities Are Read From
- [`param_value.Ar1Param`](https://statmodels7.github.io/parameters7/reference/param_value.Ar1Param.md)
  : Value of an AR(1) Parameter
- [`param_value.AutoregressiveParam`](https://statmodels7.github.io/parameters7/reference/param_value.AutoregressiveParam.md)
  : Value of an Autoregressive Parameter
- [`param_value.CompoundSymmetryParam`](https://statmodels7.github.io/parameters7/reference/param_value.CompoundSymmetryParam.md)
  : Value of a Compound Symmetry Parameter
- [`param_value.CorrelationParam`](https://statmodels7.github.io/parameters7/reference/param_value.CorrelationParam.md)
  : Value of a Correlation Parameter
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
- [`power_derivs()`](https://statmodels7.github.io/parameters7/reference/power_derivs.md)
  : Derivatives of a Power, for Composition
- [`print.parameter`](https://statmodels7.github.io/parameters7/reference/print.parameter.md)
  : Print a Constrained Parameter
- [`readable_diagonal()`](https://statmodels7.github.io/parameters7/reference/readable_diagonal.md)
  : Quantities That Are Separate Links of Separate Free Values
- [`scaled_dlog()`](https://statmodels7.github.io/parameters7/reference/scaled_dlog.md)
  : Higher Derivatives of a Scaled Log-Pseudo-Determinant
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
- [`tagged_name()`](https://statmodels7.github.io/parameters7/reference/tagged_name.md)
  : Name a Coordinate After Its Link
- [`tm_derivative()`](https://statmodels7.github.io/parameters7/reference/tm_derivative.md)
  : Derivative Components of a Transition Matrix Parameter
- [`tm_positions()`](https://statmodels7.github.io/parameters7/reference/tm_positions.md)
  : Row and Chart Coordinate of Each Free Value
- [`tuple_indices()`](https://statmodels7.github.io/parameters7/reference/tuple_indices.md)
  : The Index Tuples of a Given Width
