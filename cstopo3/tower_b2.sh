#!/bin/bash
# ---------------------------------------------------------------------------
# tower_b2.sh -- the b2 programme, sized for a 10-core / 32 GB workstation.
#
# This REPLACES slurm_array.sh for practical purposes.  The cluster array was
# written to push `cap` to 18-20 at |L| = 250, which costs ~800 s per run at
# cap 18 and still returns b2 = 1 rather than 0.  That was the wrong knob.
#
# The right knob is |L|.  At a COARSER landmark set the complex is more solid
# in dimension 3 at the same cap, because r/s is larger.  Measured:
#
#     S^3, N = 2e4, |L| = 250 :  cap 10..18 -> b2 = 42, 26, 6, 3, 1   (828 s)
#     S^3, N = 4e4, |L| = 150 :  cap 12     -> b2 = 0  EXACT, 3/3      (45 s)
#
# Whole script: well under an hour on ten cores, peak ~0.3 GB per process.
# ---------------------------------------------------------------------------
set -euo pipefail
P=9                     # leave one core for the desktop
mkdir -p results_b2/static results_b2/plateau results_b2/throat results_b2/final

echo "############ 0. self-test ############"
python3 test_all.py

echo
echo "############ 1. b2 converges for the static slices ############"
# The correct answers are S^3 -> (1,0,0) and S^1xS^2 -> (1,1,1).
# Both come out exact here; ~45-80 s per run.
python3 run_b2.py --geoms S3 S1xS2 --N 40000 --nL 150 --caps 12 \
    --seeds 8 --procs $P --out results_b2/static/

echo
echo "############ 2. the resolution plateau in |L| ############"
# Shows that b2 improves as |L| DECREASES at fixed cap, which is the point.
for nl in 120 150 200 250; do
  echo "--- |L| = $nl ---"
  python3 run_b2.py --geoms S3 --N 40000 --nL $nl --caps 12 14 \
      --seeds 4 --procs $P --out results_b2/plateau/
done

echo
echo "############ 3. the throat-radius window for the 2-cycle ############"
# The 2-cycle is the throat's S^2, and it is filled in when the throat is thin:
# lam <= 0.55 gives b2 = 0 (a solid torus, beta = (1,1,0)), lam = 0.70-0.85
# gives the correct (1,1,1), lam = 1.0 over-counts.  This is the b2 analogue of
# the short-fat-handle failure already known for b1, and it is geometry rather
# than a bug.
for lam in 0.40 0.55 0.70 0.85 1.00; do
  echo "--- lam = $lam ---"
  python3 run_cobordism.py --N 80000 --seeds 4 --nL 150 --cap 12 \
      --lam $lam --neck 1.5 --maxdim 2 --procs $P --out results_b2/throat/ | tail -3
done

echo
echo "############ 4. THE RESULT: Delta beta measured directly ############"
# Full Betti vector on both sides of the transition, no Poincare duality.
# N = 2e5 is needed because a cobordism splits the slab, so each slice sees
# only a fraction of N; at N = 8e4 the BEFORE slice is still under-resolved
# (b2 = 1..12) while the after slice is already exact.
# ~206 s per seed single-core, 0.31 GB peak.
python3 run_cobordism.py --N 200000 --seeds 32 --nL 150 --cap 12 \
    --lam 0.80 --neck 1.5 --maxdim 2 --procs $P --out results_b2/final/
python3 run_cobordism.py --merge results_b2/final/

echo
echo "############ 5. summary ############"
python3 run_b2.py --merge results_b2/static/
python3 run_b2.py --merge results_b2/plateau/
echo "--- throat sweep ---"
python3 run_cobordism.py --merge results_b2/throat/

echo
echo "Done.  Expect step 4 to report before=(1,0,0) after=(1,1,1), i.e."
echo "Delta beta = (0,+1,+1,0) measured in every component."
