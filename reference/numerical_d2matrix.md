# Numerical Second Derivatives of a Structure's Matrix

One central difference of the analytic first derivatives where a
structure supplies them, and a single mixed stencil on
[`struct_matrix`](https://statmodels7.github.io/covstructs7/reference/struct_matrix.md)
where it does not.

## Usage

``` r
numerical_d2matrix(s, eta)
```

## Arguments

- s:

  A
  [`covstruct`](https://statmodels7.github.io/covstructs7/reference/covstruct.md)
  object.

- eta:

  A numeric vector of free values.

## Value

A named list of symmetric matrices, keyed as `struct_pair_names(s)`.

## Details

Differentiating the analytic first derivative costs one
finite-difference layer rather than two, which is the rule the whole
toolkit follows: never compose differences in the same variable. When
the first derivatives are themselves numerical the two differences act
on different components for an off-diagonal pair, so they commute into a
four-point mixed stencil rather than compounding; the diagonal pairs use
the three-point second-difference stencil directly, which is again one
layer and not two.

## See also

[`struct_d2matrix`](https://statmodels7.github.io/covstructs7/reference/struct_d2matrix.md)
