#!/bin/sh
# Gate 1 main run (after controls passed): triangle U[0,1], 2000/q for q<=10, 1000/q above;
# contrast: triangle Exp(1), 500/q at every q. q=2 U[0,1] already produced by run_controls.sh.
PY=/Users/maxf/miniforge3/bin/python3
cd "$(dirname "$0")"
for q in 3 4 5 6 8 10; do $PY gate1.py --graph triangle --dist unif --q $q --n 2000; done
for q in 12 16 20; do $PY gate1.py --graph triangle --dist unif --q $q --n 1000; done
for q in 2 3 4 5 6 8 10 12 16 20; do $PY gate1.py --graph triangle --dist exp --q $q --n 500; done
