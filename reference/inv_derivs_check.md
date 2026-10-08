# The Contract Check of the Inverse's Derivatives

Refuses a parameter that is not a matrix or has no inverse, naming the
function that was called.

## Usage

``` r
inv_derivs_check(s, what)
```

## Arguments

- s:

  A parameter.

- what:

  The name of the function, for the message.

## Value

`NULL`, invisibly, or an error.
