# Log-Determinant Derivatives in the Weights

Evaluates the cyclic trace expansion \\(-1)^{n-1}\sum\_{\sigma}
\operatorname{tr}(M^{-1}P\_{\sigma(1)}\cdots M^{-1}P\_{\sigma(n)})\\ for
one tuple of weight indices, the sum running over the \\(n-1)!\\
orderings of the tail with the first index held.

## Usage

``` r
sum_struct_trace_term(minv, comp, t)
```

## Arguments

- minv:

  The inverse of the assembled matrix, computed once by the caller and
  reused for every component of the order.

- comp:

  The list of fixed components.

- t:

  An integer vector of weight indices, with repeats, of length 1 to 4.

## Value

A single number.

## Details

Holding the first index is what leaves the orderings cyclic: a cyclic
permutation leaves the trace unchanged, so fixing one position counts
each distinct cycle once. The orderings come from
[`multiset_orderings()`](https://statmodels7.github.io/parameters7/reference/multiset_orderings.md)
and are counted with multiplicity, which see.

This is a derivative in the **weights**, not in the free values;
[`sum_struct_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_logdet_derivs.md)
carries it onto the free scale.

## See also

[`sum_struct_logdet_derivs()`](https://statmodels7.github.io/parameters7/reference/sum_struct_logdet_derivs.md),
the only caller, and
[`sum_struct()`](https://statmodels7.github.io/parameters7/reference/sum_struct.md)
for the expansion.
