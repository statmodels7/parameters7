# A Short Name for a Link

The word a free name uses to say which transformation produced a
coordinate, given the link that produced it.

## Usage

``` r
link_tag(link)
```

## Arguments

- link:

  A linkfunctions7 link.

## Value

A single character string, empty for the identity link.

## Details

The tag comes from the link's class rather than from its `link_name`,
because a parametric link names itself with its parameters –
`"bounded(lwr=-0.25, upr=1)"` – and that cannot appear inside an
identifier. The identity link has no tag, so a coordinate that is
already free keeps the plain name of the quantity. A bounded link is
tagged by the transformation it performs: a doubly bounded one is a
scaled logit, and a singly bounded one a shifted logarithm. A link
written outside linkfunctions7 falls back on its own name reduced to
lowercase letters, digits and underscores.

## See also

[`tagged_name`](https://statmodels7.github.io/parameters7/reference/tagged_name.md)
