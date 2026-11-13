#!/bin/sh
# Gate 1 controls: path (C=0) U[0,1] at every q, 500 instances; triangle U[0,1] q=2 (main-run seeds, 2000).
PY=/Users/maxf/miniforge3/bin/python3
cd "$(dirname "$0")"
for q in 2 3 4 5 6 8 10 12 16 20; do $PY gate1.py --graph path --dist unif --q $q --n 500; done
$PY gate1.py --graph triangle --dist unif --q 2 --n 2000
