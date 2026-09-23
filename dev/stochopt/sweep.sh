#!/usr/bin/env bash
# Every model x seed x phase as its own process, 18 at a time.
cd ~/stoch
mkdir -p sweep1/logs
for ph in pbatch endgame; do for m in ordinal bigp panel long nonlin small; do for s in 1 2 3; do
  echo "$m $s $ph"; done; done; done |
  xargs -P 18 -L 1 bash -c 'timeout 14400 Rscript harness.R $0 $1 sweep1 $2 > sweep1/logs/$2-$0-$1.log 2>&1; echo "EXIT $?" >> sweep1/logs/$2-$0-$1.log'
echo SWEEPDONE > sweep1/DONE
