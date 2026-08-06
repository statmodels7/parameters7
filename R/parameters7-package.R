#' @keywords internal
#' @importFrom Rcpp sourceCpp
#' @useDynLib parameters7, .registration = TRUE
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' @noRd
.onLoad <- function(...) {
  S7::methods_register()
}
