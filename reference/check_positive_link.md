# Refuse a Link That Does Not Reach the Positive Half Line

Checks that an object is a linkfunctions7 link whose bounds are
contained in \\(0, \infty)\\.

## Usage

``` r
check_positive_link(link)
```

## Arguments

- link:

  The object to check.

## Value

Invisibly `TRUE`; raises an error otherwise.

## Details

A diagonal entry of a positive definite matrix is positive, so a link
onto the whole real line would let a caller build a matrix outside the
set silently. Asking the link for its declared bounds catches it at
construction, which is the only place it can be caught.
