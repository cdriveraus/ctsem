#!/usr/bin/env bash
# Every model x seed x phase as its own process, 18 at a time. Baselines first,
# so the two phases read a cached baseline instead of both computing it.
#   STOCHOPT_DIR=~/stoch/v3 ./sweep.sh <outdir>
cd "${STOCHOPT_DIR:?}"
out="${1:?outdir}"
mkdir -p "$out/logs"
run() {
  xargs -P 18 -L 1 bash -c 'timeout 14400 Rscript harness.R $0 $1 '"$out"' $2 > '"$out"'/logs/$2-$0-$1.log 2>&1; echo "EXIT $?" >> '"$out"'/logs/$2-$0-$1.log'
}
for m in ordinal bigp panel long nonlin small; do for s in 1 2 3; do
  echo "$m $s base"; done; done | run
for ph in pbatch endgame; do for m in ordinal bigp panel long nonlin small; do for s in 1 2 3; do
  echo "$m $s $ph"; done; done; done | run
echo SWEEPDONE > "$out/DONE"
