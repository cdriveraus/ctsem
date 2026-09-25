#!/usr/bin/env bash
# The dev1 half of ship.sh, copied across by it on every call. Not run by hand
# except to rebuild a label.
#
#   receive.sh init    <bundle> [donor tree]   the standing clone, once
#   receive.sh build   <label> <sha> [bundle]  a label's worktree and library
#   receive.sh harness <sha> [bundle]          the harness checkout
#
# Everything it writes is under ~/dev/ctsem-bench* and ~/dev/ctsemlib-bench-*.
set -euo pipefail
BENCH="$HOME/dev/ctsem-bench"
WTROOT="$HOME/dev/ctsem-bench-wt"
HARNESS="$HOME/dev/ctsem-bench-harness"

fetch() {
  local b="${1:-}"
  [ -n "$b" ] || return 0
  case "$b" in /*) ;; *) b="$HOME/$b" ;; esac
  git -C "$BENCH" fetch -q "$b" '+refs/heads/*:refs/remotes/ship/*'
  rm -f "$b"
}

# The inputs the Stan objects are compiled from. Two trees whose digests agree
# can share objects; ones that differ cannot, whatever their timestamps say.
stan_inputs() {
  (cd "$1" && find inst/stan inst/include -type f 2>/dev/null | LC_ALL=C sort | xargs md5sum 2>/dev/null) | md5sum | cut -d' ' -f1
}

install_to() {
  local tree="$1" lib="$2"
  mkdir -p "$lib"
  echo "installing $tree -> $lib"
  if ! (cd "$tree" && R CMD INSTALL --no-multiarch --no-byte-compile --no-docs -l "$lib" .) \
      > "$lib/install.log" 2>&1; then
    echo "INSTALL FAILED; tail of $lib/install.log:"; tail -30 "$lib/install.log"; exit 1
  fi
  # R CMD INSTALL can report success and leave nothing usable.
  [ -f "$lib/ctsem/DESCRIPTION" ] || { echo "install left no ctsem in $lib"; exit 1; }
  # Zero when the copied objects were reused; a Stan rebuild shows as several.
  local n
  n="$(grep -c '^g++ \|^gcc ' "$lib/install.log" || true)"
  echo "compiler invocations during install: $n"
}

verify_load() {
  local lib="$1"
  Rscript -e ".libPaths(c('$lib', .libPaths())); suppressMessages(library(ctsem)); cat('LOADED ctsem', as.character(packageVersion('ctsem')), 'from', find.package('ctsem'), '\n')"
}

case "${1:?usage: receive.sh init|build|harness ...}" in
  init)
    bundle="${2:?bundle}"; donor="${3:-}"
    case "$bundle" in /*) ;; *) bundle="$HOME/$bundle" ;; esac
    [ ! -e "$BENCH" ] || { echo "$BENCH exists"; exit 1; }
    git clone -q -b juliaFit "$bundle" "$BENCH"
    git -C "$BENCH" remote rename origin init-bundle
    rm -f "$bundle"
    echo "standing clone at $(git -C "$BENCH" rev-parse HEAD)"
    mkdir -p "$BENCH/src"
    # Generate the Stan sources from this tree with this machine's
    # rstantools, then look for a tree whose generated headers match byte for
    # byte and whose Stan inputs digest the same. Compare the .h, not the .cc:
    # the .cc are thin wrappers, identical across versions, and matching on
    # them is how a stale object cache went unnoticed before (CLAUDE.md).
    (cd "$BENCH" && Rscript -e 'invisible(rstantools::rstan_config())' > /dev/null)
    want="$(stan_inputs "$BENCH")"
    cands="$donor"
    [ -n "$cands" ] || cands="$(ls -dt "$HOME"/dev/ctsem-*/ 2>/dev/null | grep -v ctsem-bench || true)"
    chosen=""
    for d in $cands; do
      d="${d%/}"
      [ -f "$d/src/ctsem.so" ] || continue
      [ "$(stan_inputs "$d")" = "$want" ] || { echo "skip $d: different inst/stan or inst/include"; continue; }
      same=1
      for h in "$BENCH"/src/stanExports_*.h; do
        cmp -s "$h" "$d/src/$(basename "$h")" || { same=0; echo "skip $d: $(basename "$h") differs"; break; }
      done
      [ "$same" = 1 ] && { chosen="$d"; break; }
    done
    if [ -n "$chosen" ]; then
      echo "reusing the Linux objects of $chosen"
      # The donor's whole src/, with its own timestamps, so its objects stay
      # newer than their sources and make leaves them alone.
      cp -a "$chosen"/src/. "$BENCH/src/"
    else
      echo "no donor matches: building the Stan objects from scratch (about 20 minutes)"
    fi
    install_to "$BENCH" "$HOME/dev/ctsemlib-bench-standing"
    printf 'label=standing\nsha=%s\ndate=%s\n' "$(git -C "$BENCH" rev-parse HEAD)" "$(date -Is)" \
      > "$HOME/dev/ctsemlib-bench-standing/BENCH_BUILD"
    verify_load "$HOME/dev/ctsemlib-bench-standing"
    ;;
  build)
    label="${2:?label}"; sha="${3:?sha}"; bundle="${4:-}"
    [ -d "$BENCH/.git" ] || { echo "no standing clone at $BENCH; run ship.sh --init first"; exit 1; }
    fetch "$bundle"
    git -C "$BENCH" cat-file -e "$sha^{commit}"
    wt="$WTROOT/$label"; lib="$HOME/dev/ctsemlib-bench-$label"
    mkdir -p "$WTROOT"
    if [ -e "$wt" ]; then
      cur="$(git -C "$wt" rev-parse HEAD)"
      [ "$cur" = "$sha" ] || { echo "label $label is already $cur; use a new label"; exit 1; }
      echo "worktree $wt exists at $sha; reinstalling"
    else
      git -C "$BENCH" worktree add -q --detach "$wt" "$sha"
      mkdir -p "$wt/src"
      if [ "$(stan_inputs "$wt")" = "$(stan_inputs "$BENCH")" ]; then
        cp -a "$BENCH"/src/. "$wt/src/"
        echo "src/ copied from the standing clone"
      else
        echo "inst/stan or inst/include differ from the standing clone: Stan objects will be rebuilt"
      fi
    fi
    install_to "$wt" "$lib"
    printf 'label=%s\nsha=%s\ndate=%s\n' "$label" "$sha" "$(date -Is)" > "$lib/BENCH_BUILD"
    verify_load "$lib"
    echo "BUILD $label = $sha in $lib"
    ;;
  harness)
    sha="${2:?sha}"; bundle="${3:-}"
    fetch "$bundle"
    git -C "$BENCH" cat-file -e "$sha^{commit}"
    if [ -e "$HARNESS" ]; then
      git -C "$HARNESS" checkout -q --detach "$sha"
    else
      git -C "$BENCH" worktree add -q --detach "$HARNESS" "$sha"
    fi
    echo "HARNESS at $(git -C "$HARNESS" rev-parse HEAD): $HARNESS/dev/optimbench"
    ;;
  *) echo "unknown command $1"; exit 1 ;;
esac
