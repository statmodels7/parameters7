# Two questions that no other check asks, and that both cost a red CI push to
# learn the hard way in the sibling packages.

test_that("every object in the namespace has a help topic", {
  # The question that found seven man pages which existed on Linux and not on
  # Windows, because two topics differing only in case are one file there.
  # Asked here, on the machine where the evidence is destroyed.
  skip_if_not_installed("tools")

  db <- tryCatch(tools::Rd_db("parameters7"), error = function(e) NULL)
  skip_if(is.null(db) || !length(db), "package not installed")

  aliases <- unlist(lapply(db, function(rd) {
    tags <- vapply(rd, function(x) attr(x, "Rd_tag"), character(1))
    unlist(lapply(rd[tags == "\\alias"], function(a) trimws(paste(unlist(a),
      collapse = ""
    ))))
  }), use.names = FALSE)

  objs <- getNamespaceExports("parameters7")
  missing <- setdiff(objs, aliases)
  expect_identical(missing, character(0))
})


test_that("no two man pages differ only in case", {
  # On a case-insensitive filesystem the second overwrites the first, and
  # nothing complains: R reads the topic from \name{} inside the file, not from
  # the file name. Caught here, where both files can still exist.
  path <- system.file("man", package = "parameters7")
  if (!nzchar(path)) path <- "../../man"
  skip_if_not(dir.exists(path), "no man directory")

  files <- list.files(path, pattern = "\\.Rd$")
  skip_if(length(files) == 0L, "no Rd files")
  expect_identical(anyDuplicated(tolower(files)), 0L)
})


test_that("every exported topic is in the pkgdown reference index", {
  # pkgdown enforces this in CI and not locally, so a missing line of YAML is a
  # red push minutes later. Asking here costs a second.
  skip_if_not_installed("pkgdown")
  path <- "../.."
  skip_if_not(file.exists(file.path(path, "_pkgdown.yml")), "not the source tree")
  expect_no_error(pkgdown::check_pkgdown(path))
})


rd_info <- function(dir) {
  lapply(list.files(dir, pattern = "[.]Rd$", full.names = TRUE), function(f) {
    rd <- tools::parse_Rd(f)
    tags <- vapply(rd, function(x) attr(x, "Rd_tag"), character(1))
    grab <- function(tag) unlist(lapply(rd[tags == tag], function(x)
      trimws(paste(unlist(x), collapse = ""))))
    list(file = basename(f), name = grab("\\name"), aliases = grab("\\alias"),
         value = "\\value" %in% tags, example = "\\examples" %in% tags)
  })
}

test_that("every exported topic has a value and an executable example", {
  # The two things a first CRAN submission is most often sent back for, and
  # neither is raised by R CMD check locally. Same guard as the sibling
  # packages: a habit that lives in one package out of six is not a habit.
  dir <- if (dir.exists("../../man")) "../../man" else NULL
  skip_if(is.null(dir), "man/ not reachable from here")
  info <- rd_info(dir)
  exported <- getNamespaceExports(asNamespace("parameters7"))

  expect_equal(setdiff(exported, unlist(lapply(info, `[[`, "aliases"))),
               character())

  # the package overview page describes the package, not a return value
  pkg_page <- vapply(info, function(i) grepl("-package[.]Rd$", i$file), logical(1))
  is_exp <- vapply(info, function(i) any(i$aliases %in% exported), logical(1)) & !pkg_page
  no_value <- vapply(info, function(i) !i$value, logical(1)) & !pkg_page
  no_example <- is_exp & vapply(info, function(i) !i$example, logical(1))

  expect_equal(vapply(info[no_value], `[[`, character(1), "file"), character())
  expect_equal(vapply(info[no_example], `[[`, character(1), "file"),
               character())
})
