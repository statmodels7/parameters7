# Name a Coordinate After Its Link

Prefixes the name of a quantity with the tag of the link that carries it
onto the free scale, so `tagged_name(log_link(), "d1")` is `"log_d1"`
and `tagged_name(rhobit_link(), "rho")` is `"z_rho"`. Under the identity
link the name is left alone, the coordinate then **being** the quantity.

Every constructor in the package that takes a link builds its
`free_names` through this, which is how the convention stays uniform
across the families.

## Usage

``` r
tagged_name(link, quantity)
```

## Arguments

- link:

  A linkfunctions7 link.

- quantity:

  A character vector of quantity names, of any length.

## Value

A character vector the same length as `quantity`:
`paste0(tag, "_", quantity)` where the link has a tag, and `quantity`
unchanged where it does not.

## See also

[`link_tag()`](https://statmodels7.github.io/parameters7/reference/link_tag.md)
for the tags, and
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md),
whose `free_names` documentation states the convention and why it
matters outside the family.
