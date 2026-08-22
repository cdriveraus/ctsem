#!/bin/sh
# Refresh the vendored Julia engine in inst/julia/ContinuousTimeSEM/ from a
# ContinuousTimeSEM checkout, and record where it came from.
#
# Why vendored at all: ctJuliaSetup() used to `Pkg.add` the engine from a
# private GitLab URL pinned in inst/julia/engine.json, which meant
# backend='julia' could not be installed by anyone outside that GitLab project,
# needed network access at setup time, and pinned a *branch* rather than a
# commit. A copy inside the installed R package removes all three.
#
# The vendored tree is laid out exactly as a root-level Julia package
# (Project.toml + src/ + test/), so it is also what a standalone
# ContinuousTimeSEM repository or a registered Julia package would contain.
# That is deliberate -- see the "spinning the engine out" section of
# CPP-BACKEND.md. Once the upstream repository has the package at its root, this
# script can be replaced wholesale by:
#
#     git subtree pull --prefix=inst/julia/ContinuousTimeSEM <remote> <branch> --squash
#     git subtree push --prefix=inst/julia/ContinuousTimeSEM <remote> <branch>
#
# and edits made in either repository can flow both ways. Until then this script
# is the one-way equivalent, and it records the source commit so the two can
# always be compared.
#
# Usage: tools/sync-julia-engine.sh [path-to-ContinuousTimeSEM-repo]

set -eu

repo_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
source_repo=${1:-"$repo_root/../julia"}
package_subdir=ContinuousTimeSEM
target="$repo_root/inst/julia/ContinuousTimeSEM"

if [ ! -d "$source_repo/$package_subdir/src" ]; then
  echo "No Julia package at $source_repo/$package_subdir" >&2
  exit 1
fi

revision=$(git -C "$source_repo" rev-parse HEAD)
branch=$(git -C "$source_repo" rev-parse --abbrev-ref HEAD)
url=$(git -C "$source_repo" remote get-url origin 2>/dev/null || echo "unknown")
dirty=$(git -C "$source_repo" status --porcelain -- "$package_subdir" | head -1)
if [ -n "$dirty" ]; then
  echo "Refusing to sync: $source_repo/$package_subdir has uncommitted changes." >&2
  echo "Commit them first, so the recorded revision actually describes the copy." >&2
  exit 1
fi

# Delete the files rather than the directory tree. On a repository inside a
# synced folder (Dropbox, OneDrive) `rm -rf` on a directory intermittently fails
# with "Device or resource busy" because the syncing client holds a handle to
# it; removing the files leaves no stale content either way, and any empty
# directory left behind is immediately repopulated below.
if [ -d "$target" ]; then
  find "$target" -type f -delete
fi
mkdir -p "$target"
cp "$source_repo/$package_subdir/Project.toml" "$target/"
cp -r "$source_repo/$package_subdir/src" "$target/"
cp -r "$source_repo/$package_subdir/test" "$target/"
# The manifest is vendored too, so a user's engine environment instantiates the
# exact dependency versions this ctsem release was tested against instead of
# resolving its own.
cp "$source_repo/$package_subdir/Manifest.toml" "$target/" 2>/dev/null || true
# docs/ is not vendored: it is a Documenter site, not part of the package, and
# it would roughly double the size shipped in the R tarball.

cat > "$repo_root/inst/julia/engine.json" <<JSON
{
  "source": "vendored",
  "url": "$url",
  "branch": "$branch",
  "revision": "$revision",
  "subdir": "$package_subdir",
  "julia": "1.10"
}
JSON

echo "Vendored $package_subdir at $revision (branch $branch) into inst/julia/"
