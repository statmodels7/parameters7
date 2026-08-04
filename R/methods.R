#' @include parameter_class.R numerical_fallbacks.R
NULL


#' @title Print a Covariance Parameter
#' @name print.parameter
#'
#' @description
#' Prints the family, the shape of the matrix, the rank, the free values and
#' which of the derived quantities are numerical.
#'
#' @param x An object inheriting from class \code{\link{parameter}}.
#' @param ... Unused.
#'
#' @return Invisibly \code{x}.
#'
#' @examples
#' log_cholesky(3)
#' scaled_matrix(crossprod(diff(diag(6), differences = 2)))
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
