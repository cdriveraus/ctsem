# The parameters that vary over subjects, in the order of their random effects

The order the model declares them in: the \`indvarying\` column numbers
each varying parameter as \`.ctModelMatSetup()\` meets it in the
parameter table. This differs from matsetup's row order, which puts PARS
first for computation. Using the declared order keeps a parameter's
place among the random effects, and so the meaning of their correlation
parameters, the same whether or not it is written in PARS, and the same
on every integration method.

## Usage

``` r
.ctIndvaryingIndex(ms)
```

## Arguments

- ms:

  matsetup, as a data.frame.

## Value

integer vector of parameter numbers, one per random effect.
