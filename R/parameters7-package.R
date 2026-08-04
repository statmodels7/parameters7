#' @keywords internal
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' @noRd
.onLoad <- function(...) {
  S7::methods_register()
}
