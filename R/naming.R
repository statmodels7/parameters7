#' @include parameter_class.R
NULL

# A free name names a COORDINATE, not the quantity the coordinate produces.
# The distinction is invisible from inside a family and decisive from outside:
# a consumer flattens the free vector into scalars carrying identity links, so
# a name that reads as a constrained quantity puts a number in front of a
# reader on a scale that number is not on. An autocorrelation named "pacf1"
# and reported as 0.97 with an interval reaching past one is the shape of the
# mistake. The convention below is the one log_cholesky already follows: where
# a link carries the quantity onto the free scale, the name says which link.


#' A Short Name for a Link
#'
#' @description
#' The word a free name uses to say which transformation produced a
#' coordinate, given the link that produced it.
#'
#' @details
#' The tag comes from the link's class rather than from its
#' \code{link_name}, because a parametric link names itself with its
#' parameters -- \code{"bounded(lwr=-0.25, upr=1)"} -- and that cannot appear
#' inside an identifier. The identity link has no tag, so a coordinate that
#' is already free keeps the plain name of the quantity. A bounded link is
#' tagged by the transformation it performs: a doubly bounded one is a scaled
#' logit, and a singly bounded one a shifted logarithm. A link written
#' outside \pkg{linkfunctions7} falls back on its own name reduced to
#' lowercase letters, digits and underscores.
#'
#' @param link A \pkg{linkfunctions7} link.
#'
#' @return A single character string, empty for the identity link.
#'
#' @seealso \code{\link{tagged_name}}
#'
#' @keywords internal
link_tag <- function(link) {
  tag <- switch(attr(S7::S7_class(link), "name"),
    IdentityLink = "",
    LogLink = ,
    LowerBoundedLink = ,
    UpperBoundedLink = "log",
    LogitLink = ,
    DoublyBoundedLink = "logit",
    RhobitLink = "z",
    ProbitLink = "probit",
    ClogLogLink = "cloglog",
    LogLogLink = "loglog",
    CauchitLink = "cauchit",
    SqrtLink = "sqrt",
    InverseLink = "inv",
    InverseSqLink = "invsq",
    PowerLink = "power",
    SoftplusLink = "softplus",
    NULL
  )
  if (is.null(tag)) {
    tag <- gsub("^_|_$", "", gsub("_+", "_",
      gsub("[^a-z0-9]+", "_", tolower(link@link_name))))
    if (!nzchar(tag)) tag <- "f"
  }
  tag
}


#' Name a Coordinate After Its Link
#'
#' @description
#' Prefixes the name of a quantity with the tag of the link that carries it
#' onto the free scale, leaving the name alone when the link is the identity
#' and the coordinate therefore is the quantity.
#'
#' @param link A \pkg{linkfunctions7} link.
#' @param quantity A character vector of quantity names.
#'
#' @return A character vector the same length as \code{quantity}.
#'
#' @seealso \code{\link{link_tag}}
#'
#' @keywords internal
tagged_name <- function(link, quantity) {
  tag <- link_tag(link)
  if (nzchar(tag)) paste0(tag, "_", quantity) else quantity
}
