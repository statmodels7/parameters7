# Sines and Cosines of a Correlation Parameter's Angles

Returns, for every angle, the value and the first four derivatives **in
the free value** of both \\\sin\theta\\ and \\\cos\theta\\. These tables
are the whole derivative machinery of the family: an entry of \\L\\ is a
product of such factors, so differentiating it replaces each factor by
the derivative of the matching order, and nothing else has to be
derived.

## Usage

``` r
corr_tables(s, eta)
```

## Arguments

- s:

  A
  [`CorrelationParam()`](https://statmodels7.github.io/parameters7/reference/CorrelationParam.md)
  object, whose `param_params$link` is read.

- eta:

  A numeric vector of free values, of length `s@n_free`.

## Value

A list with two components, `sin` and `cos`, each a list of `s@n_free`
numeric vectors of length 5: the value at index 1 and the four
derivatives in the free value at indices 2 to 5.

## Details

Each angle depends on exactly one free value, so a table is indexed by
free value and no cross terms appear at this level. The chain from
\\\mathrm{d}/\mathrm{d}\theta\\ to \\\mathrm{d}/\mathrm{d}\eta\\ is
[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)'s,
applied to the link's own four derivatives, so the accuracy is the
link's and nothing is differenced.

## See also

[`compose4()`](https://statmodels7.github.io/parameters7/reference/compose4.md)
for the chain rule,
[`corr_dfactor()`](https://statmodels7.github.io/parameters7/reference/corr_dfactor.md)
which reads these tables, and
[`corr_logdet_chains()`](https://statmodels7.github.io/parameters7/reference/corr_logdet_chains.md),
which reads the `sin` half.
