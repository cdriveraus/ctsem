# Vignette runtime

The vignette `vignettes/uncertainty.qmd` was measured taking roughly 4-6
minutes of pure R/Julia execution time across its ~18 code chunks (two
backends x sampling/optimizing x priors on/off, plus an uncertainty-method
gallery and a state-explicit demo). This is acceptable for a one-off build
but is slower than the package's other qmd vignettes. If a future session
finds this vignette is slowing down routine R CMD check/CI runs, consider:
reducing n_subjects/tpoints further, trimming the uncertainty-method gallery
from 6 to 3-4 methods, or reducing iter counts on the HMC sampling chunks
(chunks named sample-stan, sample-julia, likelihood-only-sample,
likelihood-only-sample-julia, state-explicit-sample). Do not touch this
without first confirming actual CI/check timing is a real problem, since the
current runtime was deliberately tuned for clean, non-pathological results in
the normal chunks and observable pathology in the deliberately-broken
priors=FALSE+optimize=FALSE chunks.

See `.github/workflows/check-standard.yaml`: routine push/PR runs pass
`--no-vignettes` / `--no-build-vignettes` so this and the other vignettes are
not rebuilt on every push. Vignettes still ship with the package -- they are
just not re-executed by routine CI. Run that workflow manually
(`workflow_dispatch`, "Build vignettes" input) for a full pre-release check
that includes them, or locally with
`rcmdcheck::rcmdcheck(build_args = character())`.
