# Solve transform(raw) = target for raw, without knowing the transform.

The transform arrives as text using \`param\[k\]\`, so it is evaluated
rather than inverted symbolically – which means a transform this
function has never seen still works, and one it cannot bracket returns
NA rather than a guess.

## Usage

``` r
.ctDataStartInvert(text, target, bound = .ctDataStartRawBound)
```
