# The Sub-Multiset Codes of a Derivative Order

Enumerates, once per order, every sub-multiset of every index tuple the
order carries, together with the integer code each is stored under. A
sub-multiset is held as its COUNT VECTOR, one entry per free value, and
its code is that vector read as a base-five numeral.

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
each stands for.

## Details

The strings this replaces were 65 per cent of the cost of a fourth-order
component, measured: [`sort()`](https://rdrr.io/r/base/sort.html) and
[`paste()`](https://rdrr.io/r/base/paste.html) inside the key, not the
arithmetic. Counts need neither,
[`tabulate()`](https://rdrr.io/r/base/tabulate.html) being one C call,
and the code indexes a list directly. Five is a safe base because a
derivative order here is at most four, so no count can reach it.

## See also

[`ar_inv_derivative()`](https://statmodels7.github.io/parameters7/reference/ar_inv_derivative.md),
the only caller.
