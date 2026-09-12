# Duplication detectors for ctsem

Three dependency-free base-R scripts. Run from the `ctsem/` package root.

They exist because "duplicated" here usually does not mean "copy-pasted".
`ctLaplaceCorrect` and `ctParticleCorrect` are one algorithm written twice and
share almost no text; `.ctBackendAugmentedPopulation` and
`.ctBackendPopRegressionPopulation` take the same four arguments and return the
same list and read nothing alike. D1 and D2 find that shape. D3 finds
copy-paste, which is the smaller half.

```bash
cd ctsem
Rscript dev/duplication/siblings.R                                    # D1 + D2
Rscript dev/duplication/clones.R    R 60                              # D3, R
Rscript dev/duplication/clones-jl.R inst/julia/ContinuousTimeSEM/src 60   # D3, julia
```

`*-out.txt` beside each script is its output at `juliaFit` @ `a78ae672`, kept
so a later run can be diffed against it.

## What each one is for

| script | finds | blind to |
|---|---|---|
| `siblings.R` **D1 interface siblings** | functions with the same formal-argument set and no wrapper relation — two answers to one question, *however differently written* | anything whose duplicate takes different arguments |
| `siblings.R` **D2 route siblings** | every `if (pred) return(g(...))` fork — a complete inventory of where two implementations occupy one slot | forks written as `if/else` rather than early return |
| `clones.R` / `clones-jl.R` **D3** | copy-paste, after renaming | independently written duplicates (the majority) |

D1 and D2 are the ones aimed at *conceptual* duplication. D3 is the
copy-paste floor.

## Calibration notes — read before trusting a run

**D1** reports two kinds of false positive, both easy to filter and both
present in the current output:

- **S3 generic/method pairs** (`ctExtract.ctStanFit` vs `ctExtract.ctJuliaFit`)
  score 1.00 and are correct by design. Identify them by the `f.class` naming
  and the `S3method` lines in `NAMESPACE`.
- **Wrapper relations one hop away.** `.ctBackendSampleProcesses` vs
  `.ctBackendSampleEngine` scores 0.85, but the processes path reaches the
  engine via `.ctBackendSampleOneChain`, so the direct "does A call B" test
  misses it. Fixing this means a transitive call-graph check.

**D3** needs its token-diversity floor. Without it the largest reported
"clones" are `precompile_shapes.jl`'s captured data tables (579 tokens) and
long `list(a=a, b=b, ...)` constructions, which normalise to a handful of
repeated tokens and match each other trivially. The floor of 12 distinct
tokens per 60-token window dropped the Julia candidate count from 1083 to 645
and removed **only** false positives — every clone verified by hand survived
it. Raising it to 15 also removed a true positive
(`_run_chain`/`_continue_chain`), so 12 is the calibrated value, not a guess.

D3 keeps function-call names while normalising local variables and literals.
That asymmetry is deliberate: it is what stops every `for` loop matching every
other.

## Suggested CI use

D2's pair count is the useful ratchet — it may only go down. D1 and D3 are
better run by hand when starting work in an area, since both need the
judgement calls above.
