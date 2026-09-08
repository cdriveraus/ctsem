# poprank: the two measurements worth keeping

Both are contention-sensitive, so both interleave their conditions and measure
the baseline at each end of the order. Read the minima; a contamination-check
ratio away from 1.00 means the run measured the machine.

    NOT_CRAN=true Rscript dev/scratch/poprank/speed.R

| script | what it establishes |
|---|---|
| `speed.R` | gradient cost, regression form against full rank, at `k`/`m` of 20/10, 12/2 and 11/1 -- 1.56x to 1.71x, contamination checks at ratio 1.000 |
| `scale-k20.R` | the identification arithmetic at `k = 20`: unrestricted `nweak` 55 = v(v+1)/2 at every evaluation point, `poprank='auto'` 0, and `npar` 240 against 185 |

## Why the spelling matters more than the rank

An earlier prototype wrote the same reduced-rank model as loadings on auxiliary
standardised factor states, and measured *slower* than full rank at ranks 3 and
above (0.91x, 0.48x). This form is 1.56x-1.71x faster. The slower spelling had
the **smaller** augmented state (13 against 20) and **fewer** parameters (87
against 185), so neither of those is what governs the cost, and "cubic in the
augmented dimension" predicts the wrong sign twice over.

The likely reason is sparsity -- in the auxiliary spelling every cell referenced
every factor state, so the Jacobian is dense in those columns, where here each
basis cell references its own state exactly as full rank does and only the
regressed cells are expressions. That is a hypothesis: the numbers are
established, the explanation is not, and `ctsem_opcounts()` or a profile would
settle it. Do not quote the mechanism as though it were measured.
