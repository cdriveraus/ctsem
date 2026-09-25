#!/usr/bin/env bash
# Ship a commit of ctsem to dev1 and install it there as a bench build.
#
#   bash dev/optimbench/ship.sh <commit-or-branch> <label>
#   bash dev/optimbench/ship.sh --harness [<commit-or-branch>]   (default: optimbench)
#   bash dev/optimbench/ship.sh --init                           (once: the standing tree)
#   bash dev/optimbench/ship.sh --standing [<commit-or-branch>]  (move it; default juliaFit)
#
# Run it from any worktree of the ctsem repository: they share its refs, so any
# branch or commit of any job is visible from any of them.
#
# dev1 cannot reach this machine, so code travels as a git bundle holding only
# what the standing clone ~/dev/ctsem-bench does not already have. There it is
# fetched into refs/remotes/ship/*, checked out as a detached worktree at
# ~/dev/ctsem-bench-wt/<label>, given the standing clone's Linux src/ (so the
# Stan objects are not rebuilt: about 20 minutes saved, unless the label
# changed inst/stan or inst/include, in which case they are), and installed to
# ~/dev/ctsemlib-bench-<label>, with a BENCH_BUILD file naming label and sha.
# A label is bound to one sha for good: shipping a different commit under an
# existing label is refused, so a result can always be traced to its code.
#
# --harness checks the bench scripts themselves out at ~/dev/ctsem-bench-harness.
# The harness runs against any build; it is not part of the build under test.
#
# --init makes the standing clone from a full bundle of juliaFit, builds its
# src/ for Linux (reusing a donor tree's objects when its generated Stan
# headers are byte-identical, see receive.sh) and installs it to
# ~/dev/ctsemlib-bench-standing.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(git -C "$here" rev-parse --show-toplevel)"
host="${BENCH_HOST:-cd-dev1}"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

# The refs a bundle can be built from: a local branch that contains the
# commit. A bundle needs a ref, and a bare sha is refused ("empty bundle").
ref_for() {
  local rev="$1" sha="$2" r
  if git -C "$repo" show-ref --verify --quiet "refs/heads/$rev"; then
    echo "refs/heads/$rev"; return
  fi
  r="$(git -C "$repo" for-each-ref --contains "$sha" --format='%(refname)' refs/heads | head -1)"
  [ -n "$r" ] || { echo "no local branch contains $sha; commit it to a branch first" >&2; exit 1; }
  echo "$r"
}

# Bundle `ref`, excluding everything the standing clone already has, copy it
# and receive.sh across, and echo the remote bundle path ('' when dev1 already
# has the commit).
send() {
  local ref="$1" sha="$2" name="$3" have prereq=() h bundle
  ssh "$host" 'mkdir -p ~/dev/ctsem-bench-bundles'
  scp -q "$here/receive.sh" "$host:dev/ctsem-bench-bundles/receive.sh"
  if ssh "$host" "git -C ~/dev/ctsem-bench cat-file -e $sha^{commit} 2>/dev/null"; then
    echo ""; return
  fi
  have="$(ssh "$host" "git -C ~/dev/ctsem-bench for-each-ref --format='%(objectname)' 2>/dev/null" || true)"
  for h in $have; do
    if git -C "$repo" cat-file -e "$h^{commit}" 2>/dev/null &&
       git -C "$repo" merge-base --is-ancestor "$h" "$sha" 2>/dev/null; then
      prereq+=("^$h")
    fi
  done
  bundle="$tmp/optimbench-$name.bundle"
  git -C "$repo" bundle create "$bundle" "$ref" ${prereq[@]+"${prereq[@]}"} >&2
  scp -q "$bundle" "$host:dev/ctsem-bench-bundles/$name.bundle"
  rm -f "$bundle"
  echo "dev/ctsem-bench-bundles/$name.bundle"
}

case "${1:-}" in
  --init)
    ssh "$host" '[ ! -e ~/dev/ctsem-bench ]' || { echo "~/dev/ctsem-bench already exists on $host" >&2; exit 1; }
    ssh "$host" 'mkdir -p ~/dev/ctsem-bench-bundles'
    scp -q "$here/receive.sh" "$host:dev/ctsem-bench-bundles/receive.sh"
    bundle="$tmp/optimbench-init.bundle"
    git -C "$repo" bundle create "$bundle" refs/heads/juliaFit
    scp -q "$bundle" "$host:dev/ctsem-bench-bundles/init.bundle"
    rm -f "$bundle"
    ssh "$host" "bash ~/dev/ctsem-bench-bundles/receive.sh init dev/ctsem-bench-bundles/init.bundle ${BENCH_SRC_DONOR:-}"
    ;;
  --standing)
    rev="${2:-juliaFit}"
    sha="$(git -C "$repo" rev-parse --verify "$rev^{commit}")"
    ref="$(ref_for "$rev" "$sha")"
    b="$(send "$ref" "$sha" "standing-${sha:0:8}")"
    ssh "$host" "bash ~/dev/ctsem-bench-bundles/receive.sh standing $sha $b"
    ;;
  --harness)
    rev="${2:-optimbench}"
    sha="$(git -C "$repo" rev-parse --verify "$rev^{commit}")"
    ref="$(ref_for "$rev" "$sha")"
    b="$(send "$ref" "$sha" "harness-${sha:0:8}")"
    ssh "$host" "bash ~/dev/ctsem-bench-bundles/receive.sh harness $sha $b"
    ;;
  ""|-*)
    echo "usage: ship.sh <commit-or-branch> <label> | --harness [rev] | --standing [rev] | --init" >&2; exit 1 ;;
  *)
    rev="$1"; label="${2:?usage: ship.sh <commit-or-branch> <label>}"
    case "$label" in *[!A-Za-z0-9._-]*) echo "label may hold only letters, digits, . _ -" >&2; exit 1 ;; esac
    sha="$(git -C "$repo" rev-parse --verify "$rev^{commit}")"
    ref="$(ref_for "$rev" "$sha")"
    echo "shipping $rev ($sha, via $ref) as $label to $host"
    b="$(send "$ref" "$sha" "$label")"
    ssh "$host" "bash ~/dev/ctsem-bench-bundles/receive.sh build $label $sha $b"
    ;;
esac
