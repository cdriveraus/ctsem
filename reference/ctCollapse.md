# ctCollapse Easily collapse an array margin using a specified function.

ctCollapse Easily collapse an array margin using a specified function.

## Usage

``` r
ctCollapse(inarray, collapsemargin, collapsefunc, plyr = TRUE, ...)
```

## Arguments

- inarray:

  Input array of more than one dimension.

- collapsemargin:

  Integers denoting which margins to collapse.

- collapsefunc:

  function to use over the collapsing margin.

- plyr:

  Whether to use plyr.

- ...:

  additional parameters to pass to collapsefunc.

## Value

An array (or vector) with the `collapsemargin` dimensions removed,
holding `collapsefunc` applied over them.

## Examples

``` r
testarray <- array(rnorm(900,2,1),dim=c(100,3,3))
ctCollapse(testarray,1,mean)
```
