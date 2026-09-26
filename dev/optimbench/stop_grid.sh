#!/usr/bin/env bash
# Stop a grid started by run_grid.sh, and only that grid.
#
#   bash stop_grid.sh <results dir>
#
# The runner is a session leader (setsid) and everything it starts shares its
# process group, so one signal to the group stops the runner, xargs, every
# running cell and, normally, their Julias. Anything left over is found by the
# pids the cells' own harness recorded at spawn, and killed only if the
# process at that pid is still a harness or a Julia. Never by name or age:
# JuliaConnectoR's Julias are reparented to init, and on dev1 matching them by
# name or age has taken down another session's suite.
set -uo pipefail
out="${1:?usage: stop_grid.sh <results dir>}"
pidf="$out/GRID.pid"
[ -f "$pidf" ] || { echo "no GRID.pid in $out" >&2; exit 1; }
pg="$(cat "$pidf")"
if kill -0 "$pg" 2>/dev/null; then
  kill -- -"$pg" 2>/dev/null && echo "sent TERM to process group $pg"
else
  echo "runner $pg is not running"
fi
sleep 5
for f in "$out"/pids/*.R "$out"/pids/*.julia; do
  [ -f "$f" ] || continue
  p="$(cat "$f")"
  if ps -o args= -p "$p" 2>/dev/null | grep -Eq 'harness\.R|julia'; then
    kill "$p" 2>/dev/null && echo "stopped $(basename "$f") (pid $p)"
  fi
  rm -f "$f"
done
echo "stopped; finished cells keep their .rds, and run_grid.sh on the same grid resumes"
