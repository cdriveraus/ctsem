#!/usr/bin/env bash
# Run a grid of bench cells against one installed build, one R process per
# cell, at most <parallel> at a time. Written for dev1.
#
#   bash run_grid.sh <grid.txt> <label> [parallel=16]
#
# Returns at once: the runner detaches itself (setsid, stdin from /dev/null, so
# an ssh call that launches it does not hang) and writes its pid to
# <results>/GRID.pid, where
#
#   <results> = ${BENCH_RESULTS:-~/dev/ctsem-bench-results}/<label>/<grid name>
#
# The build is the library ship.sh installed, ~/dev/ctsemlib-bench-<label>
# (BENCH_LIB overrides). Data come from ~/dev/ctsem-bench-data (BENCH_DATA),
# simulated there by prepare.R on first use and read from there ever after.
#
# In <results>: one <id>.rds per finished cell, logs/<id>.log, status.tsv (one
# line per cell: id, exit code, seconds, 1-min load at the end), grid.log, and
# GRID_DONE when every cell has finished. A cell whose .rds exists is skipped,
# so launching the same grid again resumes it.
#
# Stop the whole grid with stop_grid.sh <results>, never with pkill: the
# runner, its cells and their Julias share one process group, and the cells'
# Julia pids are recorded at spawn in <results>/pids.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ "${1:-}" != "--foreground" ]; then
  grid="$(readlink -f "${1:?usage: run_grid.sh <grid.txt> <label> [parallel]}")"
  label="${2:?usage: run_grid.sh <grid.txt> <label> [parallel]}"
  par="${3:-16}"
  [ -f "$grid" ] || { echo "no such grid: $grid" >&2; exit 1; }
  name="$(basename "$grid" .txt)"
  out="${BENCH_RESULTS:-$HOME/dev/ctsem-bench-results}/$label/$name"
  mkdir -p "$out/logs" "$out/pids"
  if [ -f "$out/GRID.pid" ] && kill -0 "$(cat "$out/GRID.pid")" 2>/dev/null; then
    echo "already running: pid $(cat "$out/GRID.pid") -> $out" >&2
    exit 1
  fi
  rm -f "$out/GRID_DONE"
  setsid nohup bash "$here/run_grid.sh" --foreground "$grid" "$label" "$par" "$out" \
    >> "$out/grid.log" 2>&1 < /dev/null &
  echo $! > "$out/GRID.pid"
  sleep 2
  if kill -0 "$(cat "$out/GRID.pid")" 2>/dev/null; then
    echo "launched: pid $(cat "$out/GRID.pid") -> $out"
    echo "log: $out/grid.log"
  else
    echo "exited immediately; log:" >&2
    tail -20 "$out/grid.log" >&2
    exit 1
  fi
  exit 0
fi

shift
grid="$1"; label="$2"; par="$3"; out="$4"
# The cells run from a copy of the harness taken now. Rscript reads a script
# as it executes it, so a harness checkout moved under a running grid
# (ship.sh --harness) would corrupt every cell still going; the copy also
# records exactly which harness produced these results.
src="$here"
here="$out/harness"
mkdir -p "$here"
cp "$src"/*.R "$src"/*.jl "$src"/*.sh "$src"/*.csv "$here"/
harness_sha="$(git -C "$src" rev-parse HEAD 2>/dev/null || echo unknown)"
git -C "$src" status --porcelain -- . 2>/dev/null | grep -q . && harness_sha="$harness_sha+dirty"
echo "$harness_sha" > "$here/HARNESS_SHA"
export BENCH_DIR="$here" BENCH_LABEL="$label"
export BENCH_LIB="${BENCH_LIB:-$HOME/dev/ctsemlib-bench-$label}"
export BENCH_DATA="${BENCH_DATA:-$HOME/dev/ctsem-bench-data}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
cap="${BENCH_TIMEOUT:-7200}"
refcap="${BENCH_REFCAP:-5400}"
export BENCH_TIMEOUT="$cap" BENCH_REFCAP="$refcap"
# The harness stops its own fit at `cap` and its references at `refcap`; this
# is the backstop for a process stuck inside one Julia call, where R never gets
# to check its own limit. Two fits' worth, since the warm-up is the fit itself
# run once first (harness.R section 3).
hard=$((2 * cap + refcap + 1800))

echo "GRID $grid"
echo "harness $harness_sha (copied from $src)"
echo "label $label parallel $par results $out"
echo "build $BENCH_LIB: $(tr '\n' ' ' < "$BENCH_LIB/BENCH_BUILD" 2>/dev/null || echo 'no BENCH_BUILD')"
# nproc --all: plain nproc honours OMP_NUM_THREADS, set to 1 above.
echo "start $(date -Is) on $(hostname), $(nproc --all) cpus, load $(cat /proc/loadavg)"
echo "busiest processes at start:"
ps -eo user,pid,etime,pcpu,comm --sort=-pcpu | head -8
[ -f "$BENCH_LIB/ctsem/DESCRIPTION" ] || { echo "no build installed at $BENCH_LIB"; exit 1; }
cp "$grid" "$out/grid.txt"

# One process first: store every cell's data, and load this build's engine
# once, so the cells do not all race to precompile it.
if ! Rscript "$here/prepare.R" "$grid" > "$out/logs/_prepare.log" 2>&1; then
  echo "prepare failed:"; tail -30 "$out/logs/_prepare.log"; exit 1
fi
tail -4 "$out/logs/_prepare.log"

runcell() {
  local id="$1"; shift
  if [ -f "$out/$id.rds" ]; then echo "SKIP $id (has a result)"; return 0; fi
  local t0 code jp
  t0=$(date +%s)
  echo "START $id $(date +%H:%M:%S)"
  # No `set -e` in here: a nonzero status from this function is how xargs
  # decides to stop feeding the rest of the grid (255 stops it outright).
  timeout --foreground -k 60 "$hard" Rscript "$here/harness.R" "$id" "$@" "$out" \
    > "$out/logs/$id.log" 2>&1
  code=$?
  jp="$out/pids/$id.julia"
  if [ "$code" -eq 124 ] || [ "$code" -eq 137 ]; then
    # Killed from outside: its Julia may still be computing. The pid was
    # written by this cell's own harness when it started Julia.
    if [ -f "$jp" ] && ps -o args= -p "$(cat "$jp")" 2>/dev/null | grep -q julia; then
      kill "$(cat "$jp")" 2>/dev/null || true
    fi
  fi
  rm -f "$jp" "$out/pids/$id.R"
  printf '%s\t%s\t%s\t%s\n' "$id" "$code" "$(( $(date +%s) - t0 ))" \
    "$(cut -d' ' -f1 /proc/loadavg)" >> "$out/status.tsv"
  echo "END $id exit $code after $(( $(date +%s) - t0 ))s"
  return 0
}
export -f runcell
export out here hard

grep -v '^[[:space:]]*#' "$grid" | grep -v '^[[:space:]]*$' | \
  xargs -P "$par" -L 1 bash -c 'runcell "$@"' _
echo "GRID_DONE $(date -Is) load $(cat /proc/loadavg)"
date -Is > "$out/GRID_DONE"
