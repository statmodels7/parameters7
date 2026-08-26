# A Short Name for a Link

The word a free name uses to say which transformation produced a
coordinate, given the link that produced it.

## Usage

``` r
link_tag(link)
```

## Arguments

- link:

  A linkfunctions7 link, or any object with a `link_name` property.

## Value

A single character string, empty for the identity link. The tags in use
are `"log"` (log and either singly bounded link), `"logit"` (logit and
the doubly bounded link), `"z"` (rhobit), and `"probit"`, `"cloglog"`,
`"loglog"`, `"cauchit"`, `"sqrt"`, `"inv"`, `"invsq"`, `"power"` and
`"softplus"` for the rest. A link from outside linkfunctions7 gets its
own `link_name` reduced to lowercase letters, digits and underscores, or
`"f"` if nothing survives that.

## Details

The tag comes from the link's class, never from its `link_name`, because
a parametric link names itself with its parameters, as in
`"bounded(lwr=-0.25, upr=1)"`, and that cannot appear inside an
identifier. The identity link has no tag, so a coordinate that is
already free keeps the plain name of the quantity. A bounded link is
tagged by the transformation it performs: a doubly bounded one is a
scaled logit, and a singly bounded one a shifted logarithm. A link
written outside linkfunctions7 falls back on its own name reduced to
lowercase letters, digits and underscores.

## See also

[`tagged_name()`](https://statmodels7.github.io/parameters7/reference/tagged_name.md),
which prefixes a quantity with the tag, and
[`parameter()`](https://statmodels7.github.io/parameters7/reference/parameter.md)
for the naming convention this serves.
