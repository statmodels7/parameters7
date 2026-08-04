# Name a Coordinate After Its Link

Prefixes the name of a quantity with the tag of the link that carries it
onto the free scale, leaving the name alone when the link is the
identity and the coordinate therefore is the quantity.

## Usage

``` r
tagged_name(link, quantity)
```

## Arguments

- link:

  A linkfunctions7 link.

- quantity:

  A character vector of quantity names.

## Value

A character vector the same length as `quantity`.

## See also

[`link_tag`](https://statmodels7.github.io/parameters7/reference/link_tag.md)
