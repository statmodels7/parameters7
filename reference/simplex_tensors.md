# Derivative Tensors of the Softmax Map

The arrays \\D_1\[a, b\]\\ to \\D_4\[a, b, c, d, e\]\\ of derivatives of
\\\pi\\ in the free values, built by applying the product rule to
\\\partial_b \pi_a = \pi_a(\delta\_{ab} - \pi_b)\\ as many times as the
order asks. The first index runs over the \\K\\ categories of the value,
the rest over the \\K - 1\\ free values.

## Usage

``` r
simplex_tensors(pi_full, order = 4L)
```

## Arguments

- pi_full:

  The value, a vector of length \\K\\.

- order:

  The highest order wanted, 1 to 4.

## Value

A list with elements `d1` to `d{order}`.
