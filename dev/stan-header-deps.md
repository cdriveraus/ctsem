# Stan-generated headers and stale objects

`configure` / `configure.win` run `rstantools::rstan_config()` on every source
install. That regenerates `src/stanExports_*.{cc,h}` from `inst/stan/*.stan`
and rewrites `src/Makevars` / `src/Makevars.win` from rstantools' own fixed
template.

The staleness trap: `stanExports_X.cc` is a thin, model-name-only wrapper
around `stanExports_X.h` (an `Rcpp::exposeClass()` shim). Editing a `.stan`
file's data or parameters block regenerates `X.h` with new content and a new
mtime, but usually leaves `X.cc` byte-identical -- rstantools' own
content-diffing (`rstantools:::.add_stanfile`) then leaves `X.cc`'s mtime
untouched. `make`'s only rule for the object is R's compile rule for `.cc` ->
`.o` (an old-style *suffix* rule, `.cc.o:`, defined in
`$(R_HOME)/etc$(R_ARCH)/Makeconf` -- not a GNU *pattern* rule), so it sees no
change and relinks the stale `stanExports_X.o` against the old header. The
symptom is a data-dimension mismatch at fit time that reads exactly like a
code bug, not a build one -- there is no compile error, because nothing
recompiled.

## What it does now

`configure` and `configure.win` compare each generated header against its
object, right after `rstan_config()` has regenerated the headers, and delete
the object when the header is newer:

```sh
if [ -n "$(find "$_ctsem_hdr" -newer "$_ctsem_obj" 2>/dev/null)" ]; then
  rm -f "$_ctsem_obj"
fi
```

`make` then rebuilds it because it is missing. The model list comes from
`inst/stan/*.stan` at configure time, so an added or removed `.stan` file
needs no matching edit. When it fires it says so on stdout, which appears in
the install log.

The important property is what it does *not* do: it leaves the generated
`src/Makevars`/`Makevars.win` byte-identical to rstantools' template. That
keeps `rstan_config()` from rewriting them (it writes only on a content
difference), so their mtime never moves, and it keeps pkgbuild's own staleness
check meaningful rather than permanently tripped -- pkgbuild still precleans
when a real header, such as one in `inst/include`, is newer than the built
dll.

## Do not move this back into Makevars

Three versions tried that, and each failed differently. All three were
reproduced directly rather than reasoned about, and all three are silent.

1. **A `%` pattern rule does nothing.** `stanExports_%.o: stanExports_%.h`
   looks equivalent and is the obvious thing to write. A prerequisite-only GNU
   *pattern* rule is not consulted for a target that make is actually building
   via an old-style *suffix* rule (`.cc.o:`, what Makeconf defines) -- the two
   implicit-rule mechanisms don't combine for that purpose. `make` reports
   "Nothing to be done" even though the header is newer. Reproduced against
   Rtools' GNU Make 4.4.1 with a minimal two-rule Makefile.

2. **An explicit-name rule becomes make's default goal.** `tools:::.shlib_internal`
   builds `makefiles <- c("Makevars", makefiles)` and invokes make without
   naming a goal, so the FIRST explicit target in Makevars becomes the default.
   make built that single object, printed "Nothing to be done", exited 0, and
   `R CMD INSTALL` ran happily through every remaining stage before failing at
   the load check with no shared object and no error anywhere in 36,500 lines
   of log. Only a from-scratch build exposes it: a tree that already has every
   object still links after building one. Three from-scratch failures on Linux
   against a Windows tree that passed.

3. **Any edit at all makes the file get rewritten every run, and that wipes
   every object.** `rstantools:::.add_stanfile()` writes when the content
   differs from its template, which an appended rule guarantees -- same bytes
   afterwards, new mtime, every configure run. `pkgbuild` counts
   `src/Makevars*` as a *header* (`pkgbuild:::headers`), and
   `pkgbuild:::needs_clean()` is "any header newer than `src/<pkg>.dll`". When
   true, `compile_dll()` passes `--preclean`, which runs `shlib-clean`
   (`rm -f $(OBJECTS)` -- objects only, the dll is left behind) and every stan
   object is rebuilt from scratch. That is the path devtools, roxygen and
   RStudio's build all take. It is self-sustaining, because an install whose
   `make` has nothing to do does not relink the dll, leaving the Makevars that
   the same install's configure step rewrote newer than it. Observed as a full
   stan recompile on roughly every second RStudio build, with `src/` holding a
   fresh `RcppExports.o`, no `stanExports_*.o` at all, and a `ctsem.dll`
   several days old.

A fourth version kept the rules and restored the Makevars mtime with
`cp -p` / `cmp -s` / `touch -r` around the `rstan_config()` call. It worked,
but it was a third workaround propping up the first, and the check moved out
of Makevars entirely instead.

## To prove it still works after touching this

1. Run `sh configure.win` (with `R_HOME` set) two or three times on a built
   tree. `src/Makevars.win`'s mtime must not move at all, no object may be
   removed, and `pkgbuild:::needs_clean(".")` must stay `FALSE`.
2. `touch src/stanExports_ctsmgen.h`, run it again: that one object is removed
   and reported, and `stanExports_ctsm.o` is not.
3. From a tree with no `src/` directory at all, confirm the script exits 0 and
   leaves the Makevars byte-identical to rstantools' template.
4. End to end, edit an `inst/stan/*.stan` file with a **genuine code change**,
   not a `touch` or a comment. `stanc` strips comments, and rstantools'
   content-diffing means a change that doesn't alter the compiled output
   doesn't even bump the generated header's mtime -- there is nothing to
   detect a rebuild from. A harmless forced change: add an unused local, e.g.
   `transformed data{ int probe = 1; }`, to `ctsmgen.stan`, the cheaper of the
   two models. Reinstall: only that model's object and the final link should
   redo. Revert and reinstall again as a round trip.
