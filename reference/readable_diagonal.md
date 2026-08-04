# Quantities That Are Separate Links of Separate Free Values

Assembles the declaration of a family whose quantities are one scalar
link each of one free value each, so that the Jacobian is the diagonal
matrix of the inverse links' first derivatives.

## Usage

``` r
readable_diagonal(links, eta, nm, transform, label)
```

## Arguments

- links:

  A list of linkfunctions7 links, one per quantity, in the order of the
  free values they read.

- eta:

  A numeric vector of free values.

- nm:

  A character vector naming the quantities.

- transform:

  A character vector naming the scale each interval is built on.

- label:

  A single string naming the group.

## Value

A list as described in
[`param_readable`](https://statmodels7.github.io/parameters7/reference/param_readable.md).
