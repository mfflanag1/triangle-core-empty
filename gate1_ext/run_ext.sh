#!/bin/sh
# Post-gate extension run (2026-11-11): core-existence share only, triangle U[0,1].
# q=24: n=500 with MILP gap (first launch). Under machine load the q=40 MILPs did not finish
# a single instance in 5 minutes, so q=28, 32, 40 run in --lp-only mode: core nonempty iff the
# dual-simplex LP vertex is integral (a.s.-unique optimum; gate 1 saw 0 fractional-vertex/zero-gap
# cases in 21,000 triangle instances). ext.py resumes from existing output.
PY=/Users/maxf/miniforge3/bin/python3
cd "$(dirname "$0")"
$PY ext.py --q 24 --n 500 --workers 6
$PY ext.py --q 40 --n 500 --workers 6 --lp-only
$PY ext.py --q 32 --n 500 --workers 6 --lp-only
$PY ext.py --q 28 --n 500 --workers 6 --lp-only
# Second extension (same session): more q=40 and a q=48 point, LP-only, to test for a floor.
$PY ext.py --q 40 --n 2000 --workers 6 --lp-only
$PY ext.py --q 48 --n 1000 --workers 6 --lp-only
