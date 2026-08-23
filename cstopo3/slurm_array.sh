#!/bin/bash
#SBATCH --job-name=cstopo3-b2
#SBATCH --array=0-31
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/b2_%A_%a.out
# ---------------------------------------------------------------------------
# THE CLUSTER JOB: convergence of b_2.
#
# b_0 and b_1 are cheap enough for the tower (see tower_run.sh) -- do not spend
# cluster time on them.  b_2 is the one that needs it: it requires the
# 3-simplices, C(cap,4) per facet, and converges only at cap ~ 18-20.
#
# 2 geometries x 4 caps x 8 seeds = 64 runs over 32 array tasks = 2 runs each.
# Measured single-core cost per run (S^3, N = 20000, |L| = 250):
#     cap 14 ->   67 s        cap 18 -> ~550 s   (extrapolated)
#     cap 16 ->  198 s        cap 20 -> ~1400 s  (extrapolated)
# so a task holding two cap-20 runs needs about an hour.  4 h of wall is ample.
#
# Memory stays under ~2 GB at cap 18.  The --guard flag skips any run whose raw
# 3-simplex count would exceed the limit rather than dying on the node.
# ---------------------------------------------------------------------------
set -euo pipefail
mkdir -p logs results_b2

module load python/3.11 2>/dev/null || true
python3 -c "import numpy; print('numpy', numpy.__version__)"   # only dependency

srun python3 run_b2.py \
    --geoms S3 S1xS2 \
    --N 20000 --nL 250 \
    --caps 14 16 18 20 \
    --seeds 8 \
    --guard 40000000 \
    --task "${SLURM_ARRAY_TASK_ID}" --ntasks 32 \
    --out results_b2/

# When the whole array has finished:
#     python3 run_b2.py --merge results_b2/
#
# What to look for: b2 -> 0 for S3 and b2 -> 1 for S1xS2 as cap rises.  If the
# trend flattens above zero, the cover is not becoming solid in dimension 3 and
# the next knob is a SMALLER |L| at the same cap (larger balls), not a larger N
# -- see the |L| row in HANDOFF.md, where |L| = 150 at cap 12 already gives the
# correct b2 for both geometries with exact metrics.
