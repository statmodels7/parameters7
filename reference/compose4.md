# Compose Two Scalar Maps, to Fourth Order

The derivatives of \\f(g(x))\\ of orders one to four, from the
derivatives of \\f\\ at \\g(x)\\ and of \\g\\ at \\x\\.

## Usage

``` r
compose4(fd, gd)
```

## Arguments

- fd:

  A list or vector of the four derivatives of the outer map at \\g(x)\\,
  in order. Read with `[[`, so a list, a numeric vector or a length-4
  list of numeric vectors all serve. The **value** of \\f\\ is not
  needed and is not read.

- gd:

  A list or vector of the four derivatives of the inner map at \\x\\, in
  order, read the same way. The value of \\g\\ is not needed here
  either; the caller has already evaluated \\f\\'s derivatives at it.

## Value

A list of four elements, the composite derivatives in order, each the
shape the arithmetic on `fd` and `gd` produces: a single number where
both are scalar, a vector where either is, a matrix where either is.

## Details

Faa di Bruno's formula, whose coefficients are the numbers of set
partitions of a given shape: \$\$(f \circ g)' = f' g',\$\$ \$\$(f \circ
g)'' = f'' g'^2 + f' g'',\$\$ \$\$(f \circ g)''' = f''' g'^3 + 3 f'' g'
g'' + f' g''',\$\$ \$\$(f \circ g)'''' = f'''' g'^4 + 6 f''' g'^2 g'' +
3 f'' g''^2 + 4 f'' g' g''' + f' g''''.\$\$ The four coefficients of the
last line count the partitions of four elements into four singletons, a
pair and two singletons, two pairs, a triple and a singleton, and one
block.

Every argument may be a vector, in which case the composition is applied
elementwise and the result has the same length. That is how the families
use it: one call composes a whole table of angles or lags at once.

Verified against `numDeriv` on \\\exp(\sin x)\\ at \\x = 0.7\\, the
first two orders agreeing to eight figures.

## See also

[`leibniz_gram()`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md)
for the other piece of shared arithmetic,
[`power_derivs()`](https://statmodels7.github.io/parameters7/reference/power_derivs.md)
for the commonest outer map here, and
[`numericals7::set_partitions()`](https://statmodels7.github.io/numericals7/reference/set_partitions.html),
which enumerates the partitions the coefficients count.
