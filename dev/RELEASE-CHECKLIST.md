# Release checklist

Written 2026-09-03, from the things that actually went wrong during the review
rather than from a generic list. Each item exists because skipping it cost time
or shipped something wrong.

## Before anything else

**Push, and check what is already pushed.** `git rev-list --count
origin/juliaFit..juliaFit`. During the review this drifted to twenty-two
commits, and a stale remote makes every later question harder.

## The build

**Objects rebuild by themselves now; do not delete them.** Until 2026-09-03
the makefile did not declare the generated headers as prerequisites of their
objects, so a months-old object was silently relinked when a `.stan` file
changed, and the working advice was to delete `src/stanExports_*.o` before
installing. `configure` now emits an explicit rule per model, so make
rebuilds exactly the models whose headers changed and nothing else. The
deletion step is obsolete and should not be reinstated.

If you ever change how `src/` is produced, re-prove this from NOTHING rather
than from a tree that already has objects. The first version of that fix
passed a targeted rebuild on a fully built tree and broke installation from
scratch on every platform, because its rule became make's default goal: make
built one object, said "Nothing to be done", exited zero, and the install
failed only at the load step with no error anywhere in the log. A build fix
verified on a tree that already has everything built is not verified.

**Two things legitimately compile from scratch and are not faults.** `R CMD
build` runs `* cleaning src` before making a tarball, so testing a real build
costs a full recompile; keep that rare. And a fresh ship to dev1 excludes
objects and libraries deliberately, because Windows objects will not link on
Linux.

**Install into a library directory that already exists.** `R CMD INSTALL -l
../somewhere` fails immediately if the directory is absent, in a way that looks
like it ran.

## Shipped data

**Regenerate `ctstantestfit` whenever the fit object's structure changes.** It
stores `standata`, so a structural change makes anything feeding it back to stan
fail. This happened when the `stateref` column was added: the shipped fit still
carried a nine-column parameter setup array against code declaring ten, and
`ctGenerateFromFit(ctstantestfit)` errored.

**Do not regenerate `ctstantestdat` at the same time without meaning to.**
`ctdataupdate()` rebuilds both. The data is produced under a fixed seed, but
generation behaviour changes, so regenerating it moves values that many tests
compare against. Refresh the fit alone unless the data change is the point.

## Tests

**`NOT_CRAN=true` or the julia tests silently skip** and a file reporting zero
assertions looks identical to a file that passed. Check assertion counts, not
just the absence of failures.

**Do not trust a green run from `testall()` alone** for anything before
September 2026: its workers loaded whichever ctsem was installed rather than the
tree in hand, so it may have been testing something else entirely. Fixed now,
but historical results from it prove less than they appear to.

**Run the stan-julia parity suite deliberately.** It skips unless
`CTSEM_JULIA_PROJECT` is set, so it does not run in an ordinary suite. It is the
only mechanical check that the two backends agree, and several divergences this
review found by hand would have been caught by it.

**Run the julia engine suite on both platforms.** Its assertions have not
historically been portable: one test passed on Windows and failed on Linux for
months because a fixture sat where a transform's derivative had collapsed. Both
platforms now give the same count, and keeping it that way needs both to be run.

## Numbers and claims

**A fit reporting `converged = FALSE` is a finding, not noise.** Until September
2026 the saturation guard tested the magnitude of every raw parameter against a
single constant, so ordinary fits on data with a large mean reported themselves
as failures. It now tests whether a transform has stopped responding.

**Benchmark on dev1, not locally, and name the machine on every number.**
Contention alone has produced a 140% spread on identical work here.

## NEWS

Written last, in one deliberate pass, from the user's side. A bug introduced and
fixed within the same unreleased cycle gets no entry. See the house rules in
`CLAUDE.md`; the register is one line per change, no mechanism and no
measurements.
