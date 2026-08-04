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
  in order.

- gd:

  A list or vector of the four derivatives of the inner map at \\x\\, in
  order.

## Value

A list of four elements, the composite derivatives in order.

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
elementwise and the result has the same length.

## See also

[`leibniz_gram`](https://statmodels7.github.io/parameters7/reference/leibniz_gram.md)
