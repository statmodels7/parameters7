# Third Derivatives of a Correlation Parameter

Closed form, the same Leibniz rule on \\R = LL^\top\\ with the three
differentiations distributed over the two factors. Each \\\partial^S L\\
comes from
[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md),
which returns `NULL` whenever \\S\\ spans two rows of \\L\\ or reaches
past the columns an entry involves, so most terms of the sum are skipped
instead of computed and discarded.

The angles reach the free scale through a bounded link, so the chain to
third order is
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)'s
and the accuracy is the link's; nothing is differenced. The diagonal is
exactly zero.

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object.

- eta:

  A numeric vector of free values, of length `s@n_free`, already checked
  by the generic.

- ...:

  Unused, and accepted so the signature matches the generic's.

## Value

A list of `choose(s@n_free + 2, 3)` symmetric matrices keyed as
`param_tuple_names(s, 3)` and in that order, each with a zero diagonal.

## See also

[`corr_derivative()`](https://statmodels7.github.io/parameters7/reference/corr_derivative.md),
which assembles it, and
[`param_d4.CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/param_d4.CorrelationParam.md)
for the order above.
