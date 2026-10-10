# The Sub-Multiset Codes of a Derivative Order

Enumerates, once per order, every sub-multiset of every index tuple that
the order carries, together with the integer code under which each is
stored. A sub-multiset is held as its count vector, one entry per free
value, and its code is that vector read as a base-five numeral.

## Usage

``` r
ar_inv_codes(idx, n)
```

## Arguments

- idx:

  The index tuples of the order, as
  [`param_tuple_indices()`](https://statmodels7.github.io/parameters7/reference/param_tuple_indices.md)
  returns them.

- n:

  The number of free values.

## Value

A list with `pow`, the base-five weights; `codes`, one integer vector
per tuple holding the code of `t[take]` for every subset of positions in
the order the bits enumerate them; `comp`, the matching codes of the
complements; and `all`, every distinct code with `cnt`, the count vector
that each stands for.

## Details

A count vector needs neither a string key nor a sort:
[`tabulate()`](https://rdrr.io/r/base/tabulate.html) builds it in one
call, and the code indexes a list directly. Five is a safe base because
a derivative order here is at most four, so no count can reach it.

## See also

[`ar_inv_structure()`](https://statmodels7.github.io/parameters7/reference/ar_inv_structure.md),
the only caller.
