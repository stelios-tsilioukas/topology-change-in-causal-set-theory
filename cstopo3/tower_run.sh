#!/bin/bash
# ---------------------------------------------------------------------------
# NOA tower, 10 cores / 32 GB.  Everything here fits comfortably; the whole
# script is roughly 1-2 hours wall.  Run it before touching the cluster.
#
# Memory is a non-issue for b0/b1: the largest run measured used 0.10 GB.
# ---------------------------------------------------------------------------
set -euo pipefail
P=9                     # leave one core for the desktop
mkdir -p results

echo "############ 0. self-test (should be under two minutes) ############"
python3 test_all.py

echo
echo "############ 1. cost scaling on this hardware ############"
python3 run_static.py --scaling

echo
echo "############ 2. static-slice plateau, 32 seeds ############"
# The plateau in the shadow-cardinality resolution.  cap = 8 under-resolves,
# cap >= 10 is the plateau.  This is the figure for the paper.
python3 run_static.py --N 20000 --seeds 32 --nL 250 \
    --caps 6 8 10 12 14 16 --maxdim 1 --procs $P --out results/

echo
echo "############ 3. the headline: cobordism, 128 seeds ############"
# S^3 -> S^1 x S^2 detected from the causal order alone.  ~8 s per seed on one
# core, so 128 seeds across 9 processes is about two minutes.
python3 run_cobordism.py --N 20000 --seeds 128 --nL 250 --cap 14 \
    --procs $P --out results/
python3 run_cobordism.py --merge results/

echo
echo "############ 4. cobordism at larger N, 32 seeds ############"
# Confirms the detection does not degrade with resolution.  N = 160000 costs
# ~63 s per seed single-core.
python3 run_cobordism.py --N 80000  --seeds 32 --nL 350 --cap 14 --procs $P --out results/
python3 run_cobordism.py --N 160000 --seeds 32 --nL 400 --cap 14 --procs $P --out results/

echo
echo "############ 5. throat-radius and neck sweep ############"
# Maps the region of handle geometry in which the 1-cycle survives.  A short fat
# handle (large lam, small neck) has a small-feature 1-cycle and is filled in by
# the cover: that failure is geometric, not a bug, and worth reporting.
for lam in 0.25 0.40 0.55; do
  for neck in 0.0 0.75 1.50; do
    echo "--- lam=$lam neck=$neck ---"
    python3 run_cobordism.py --N 20000 --seeds 16 --lam $lam --neck $neck \
        --procs $P --out results/ | tail -2
  done
done

echo
echo "############ 6. b2 -- as far as the tower will take it ############"
# Full Betti vector including b2.  Costs minutes per run rather than seconds.
# The remainder (cap 18, 20) belongs on the cluster: see slurm_array.sh.
python3 run_b2.py --geoms S3 S1xS2 --N 20000 --nL 250 \
    --caps 10 12 14 16 --seeds 4 --procs $P --out results_b2/
python3 run_b2.py --merge results_b2/

echo
echo "Done.  results/ and results_b2/ hold the JSON; the --merge calls above"
echo "print the summaries.  Next step is slurm_array.sh for cap 18-20."
