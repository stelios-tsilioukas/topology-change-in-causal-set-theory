#!/usr/bin/env python3
"""
run_trousers.py -- end-to-end demonstration of topology change in a 2+1D
causal set, using ONLY the causal order.

Spatial slices go  one circle -> two circles, the 2+1D analogue of
S^3 -> S^1 x S^2.  The pipeline is:

    sprinkle -> causal matrix -> inextendible antichain -> MRS nerve
             -> GF(2) homology -> Betti plateau -> Morse count

Job-array parallel: pass --task K --ntasks M and each task handles the seeds
K, K+M, K+2M, ...  Results are written as one JSON per task; merge with
--merge.

Examples
--------
    # single run
    python run_trousers.py --N 8000 --seeds 8

    # SLURM array of 32 tasks, 4 seeds each
    #SBATCH --array=0-31
    python run_trousers.py --N 40000 --seeds 128 \
        --task $SLURM_ARRAY_TASK_ID --ntasks 32 --out results/

    # merge
    python run_trousers.py --merge results/
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cstopo.geometry import trousers_2p1, cylinder_2p1
from cstopo.topology import antichain_from_band, stability_scan, morse_count

DEFAULT_CAPS = (4, 6, 8, 10, 12, 14)


def one_seed(N, seed, caps, nslices, band, T=2.0, t_split=1.0, verbose=False):
    """Run the full pipeline for one sprinkling.  Returns a dict."""
    t0 = time.time()
    t, x, R, meta = trousers_2p1(N, T=T, t_split=t_split, seed=seed)
    t_build = time.time() - t0

    # foliate
    times = np.linspace(0.35 * T, 0.90 * T, nslices)
    slices = []
    for tc in times:
        A = antichain_from_band(R, t, tc, band)
        if len(A) < 8:
            slices.append(dict(t=float(tc), nA=int(len(A)), scan=None))
            continue
        scan = stability_scan(R, A, caps=caps, maxdim=2)
        slices.append(dict(
            t=float(tc), nA=int(len(A)),
            scan={int(c): (list(v[0]) if v[0] != "too_big" else "too_big")
                  for c, v in scan.items()}))

    # plateau Betti at each slice: the value taken over the most caps
    seq = []
    for s in slices:
        if s["scan"] is None:
            seq.append(None); continue
        vals = [tuple(v[:2]) for v in s["scan"].values() if v != "too_big"]
        seq.append(max(set(vals), key=vals.count) if vals else None)

    chi, events = morse_count([v for v in seq if v is not None])
    return dict(N=N, seed=int(seed), t_build=t_build,
                t_total=time.time() - t0, times=[float(v) for v in times],
                slices=slices, plateau=[list(v) if v else None for v in seq],
                morse_chi=int(chi), events=[[int(a), int(b)] for a, b in events])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=8000)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--seed0", type=int, default=1000)
    ap.add_argument("--nslices", type=int, default=9)
    ap.add_argument("--band", type=float, default=0.08)
    ap.add_argument("--caps", type=int, nargs="+", default=list(DEFAULT_CAPS))
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--merge", type=str, default=None)
    ap.add_argument("--mpi", action="store_true",
                    help="distribute seeds with mpi4py instead of a job array")
    args = ap.parse_args()

    if args.merge:
        recs = []
        for fn in sorted(os.listdir(args.merge)):
            if fn.endswith(".json"):
                recs += json.load(open(os.path.join(args.merge, fn)))
        summarise(recs)
        return

    rank, size = args.task, args.ntasks
    if args.mpi:                                    # ---- MPI hook ----
        from mpi4py import MPI                      # noqa: F401
        comm = MPI.COMM_WORLD
        rank, size = comm.Get_rank(), comm.Get_size()

    seeds = [args.seed0 + s for s in range(args.seeds)][rank::size]
    caps = tuple(args.caps)
    recs = []
    for s in seeds:
        r = one_seed(args.N, s, caps, args.nslices, args.band)
        recs.append(r)
        print(f"[rank {rank}] seed {s}: chi={r['morse_chi']:+d} "
              f"({r['t_total']:.1f}s)", flush=True)

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        with open(os.path.join(args.out, f"task{rank:04d}.json"), "w") as f:
            json.dump(recs, f)
    else:
        summarise(recs)


def summarise(recs):
    from collections import Counter
    print("\n" + "=" * 74)
    print(f" SUMMARY over {len(recs)} sprinklings")
    print("=" * 74)
    if not recs:
        return
    N = recs[0]["N"]
    early = Counter(); late = Counter(); chis = Counter()
    for r in recs:
        pl = r["plateau"]
        k = len(pl) // 2
        for v in pl[:k]:
            if v: early[tuple(v)] += 1
        for v in pl[k:]:
            if v: late[tuple(v)] += 1
        chis[r["morse_chi"]] += 1
    print(f" N = {N}")
    print(f" early slices (expect (1,1)): {dict(early.most_common(4))}")
    print(f" late  slices (expect (2,2)): {dict(late.most_common(4))}")
    print(f" Morse chi                  : {dict(chis.most_common(5))}")
    tt = np.mean([r["t_total"] for r in recs])
    print(f" mean wall time per sprinkling: {tt:.1f} s")


if __name__ == "__main__":
    main()
