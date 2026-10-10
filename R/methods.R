#' @include parameter_class.R numerical_fallbacks.R
NULL


#' @title Print a Constrained Parameter
#' @name print.parameter
#'
#' @description
#' Prints a one-screen summary of a parametrization: the family name, the shape
#' and rank of the matrix, how many free values there are and what they are
#' called, and which derived quantities (nine for a matrix family, four for
#' [simplex()] and [transition_matrix()]) come from the base class instead of a
#' closed form.
#'
#' @details
#' The matrix block is printed only for a [matrix_parameter()], because a
#' [simplex()] has no dimension or rank to report. The free names are truncated at twelve, with a count of the rest, so a
#' large unstructured covariance does not fill the screen. The last line,
#' `From the base class`, reports [param_is_numerical()]: `none` means that
#' every quantity has a method of the family's own, and a list of generics names
#' the ones that the base class supplies. Every family in this package prints
#' `none`.
#'
#' A rank-deficient family prints the dimension of its null space beside the
#' rank, which shows at once that a penalty is improper.
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
