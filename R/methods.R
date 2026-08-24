#' @include parameter_class.R numerical_fallbacks.R
NULL


#' @title Print a Constrained Parameter
#' @name print.parameter
#'
#' @description
#' Prints a one-screen summary of a parametrization: the family name, the shape
#' and rank of the matrix, the role label, how many free values there are and
#' what they are called, and which of the nine derived quantities come from the
#' base class instead of a closed form.
#'
#' @details
#' Three things about the output are worth knowing. The matrix block is printed
#' only for a [matrix_parameter()], a [simplex()] having no dimension or rank to
#' report. The free names are truncated at twelve, with a count of the rest, so a
#' large unstructured covariance does not fill the screen. And the last line,
#' `From the base class`, is [param_is_numerical()]'s answer: `none` means every
#' quantity is exact, and a list of generics names the ones a stencil computes.
#' Every family this package ships prints `none`.
#'
#' A rank-deficient family prints its null space's dimension beside the rank,
#' which is the quickest way to see that a penalty is improper.
#'
#' @param x An object inheriting from class [parameter()].
#' @param ... Unused, and accepted so the signature matches the generic's.
#'
#' @return Invisibly `x`, so a print inside a pipe does not break it.
#'
#' @seealso [param_is_numerical()], whose answer the last line reports, and
#'   [check_parameter()] for a validation report instead of a description.
#'
#' @examples
#' # A full-rank family: nothing from the base class.
#' log_cholesky(3)
#'
#' # A rank-deficient one prints its null space beside the rank.
#' scaled_matrix(crossprod(diff(diag(6), differences = 2)))
#'
#' # A family that is not a matrix prints no matrix block.
#' simplex(4)
#'
#' # Twelve free names at most, with the rest counted.
#' log_cholesky(6)
#'
#' @keywords internal
S7::method(print, parameter) <- function(x, ...) {
  cat(sprintf("Parameter: %s\n", x@param_name))
  if (S7::S7_inherits(x, matrix_parameter)) {
    cat(sprintf("Matrix:    %d x %d, symmetric\n", x@dimension, x@dimension))
    cat(sprintf(
      "Rank:      %d of %d%s\n", x@rank, x@dimension,
      if (x@rank < x@dimension) {
        sprintf(" (null space of dimension %d)", x@dimension - x@rank)
      } else {
        ""
      }
    ))
    cat(sprintf("Role:      %s\n", x@role))
  }

  cat(sprintf("\nFree values: %d\n", x@n_free))
  if (x@n_free) {
    shown <- utils::head(x@free_names, 12L)
    cat("  ", paste(shown, collapse = ", "),
      if (x@n_free > 12L) sprintf(", ... (%d more)", x@n_free - 12L) else "",
      "\n",
      sep = ""
    )
  }

  num <- param_is_numerical(x)
  cat(sprintf(
    "\nFrom the base class: %s\n",
    if (any(num)) paste(names(num)[num], collapse = ", ") else "none"
  ))

  invisible(x)
}
