#!/usr/bin/env python3
"""
run_wormhole.py -- SUPERSEDED, kept so the failure can be reproduced.

This is the first 3+1D attempt.  It carries the 2+1D shadow-cardinality
resolution (cap ~ 10) over to a three-dimensional slice unchanged, which
under-resolves it: the modal Betti vector is noise at N <= 3e4.  The diagnosis
and the fix are in ../cstopo3/ ; use that for the 3+1D case.

run_wormhole.py -- test whether a sprinkled Giddings-Strominger wormhole yields
the topology-change signature  delta_beta = (0,+1,+1,0),  hence delta_chi = -2.

The spatial slices of the GS geometry are three-spheres away from the throat.
Because chi(Sigma^3) == 0 identically, the Euler characteristic of a slice
carries no information; the discriminating observable is the Betti VECTOR,

    S^3          : (1,0,0,1)
    S^1 x S^2    : (1,1,1,1)

so that creation is  Delta beta = (0,+1,+1,0)  (Morse index 1) and destruction
its dual (index 3), giving  delta_chi = (-1)^1 + (-1)^3 = -2  by Eq. (morse).

CAUTION, please read before interpreting a run:
  * beta = (1,1,1,1) is NECESSARY but NOT SUFFICIENT for S^1 x S^2 -- the wedge
    S^1 v S^2 v S^3 has the same Betti vector.  Use --pseudomanifold to apply
    the cheap necessary check (every 2-face in exactly two 3-simplices).
  * a NEGATIVE result at accessible N is a statement about the observable, not
    about the physics.  Report the scaling, not just the outcome.

Recommended staging: run --scaling first at N = 1e4, 3e4, 1e5 and inspect how
the plateau width and detection rate move BEFORE requesting a large allocation.

Examples
--------
    python run_wormhole.py --scaling
    python run_wormhole.py --N 100000 --seeds 16 --task $SLURM_ARRAY_TASK_ID \
        --ntasks 16 --out results_wh/
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cstopo.geometry import gs_wormhole
from cstopo.topology import (antichain_from_band, stability_scan, morse_count,
                             layer_index, mrs_nerve_fast, betti_gf2)

TARGET_S3 = (1, 0, 0, 1)
TARGET_WH = (1, 1, 1, 1)


def is_pseudomanifold(maximal, d=3):
    """Necessary condition for a d-manifold: every (d-1)-face lies in exactly
    two d-simplices.  Filters the wedge S^1 v S^2 v S^3."""
    from itertools import combinations
    from collections import Counter
    top = [s for s in maximal if len(s) == d + 1]
    if not top:
        return False
    cnt = Counter()
    for s in top:
        for f in combinations(sorted(s), d):
            cnt[f] += 1
    return all(v == 2 for v in cnt.values())


def one_seed(N, seed, caps, nslices, band, a0=1.0, tau_max=2.0,
             check_pm=False):
    t0 = time.time()
    tau, x, R, meta = gs_wormhole(N, a0=a0, tau_max=tau_max, seed=seed)
    t_build = time.time() - t0

    lo, hi = tau.min(), tau.max()
    times = np.linspace(0.75 * lo, 0.75 * hi, nslices)
    out = []
    for tc in times:
        A = antichain_from_band(R, tau, tc, band)
        if len(A) < 10:
            out.append(dict(tau=float(tc), nA=int(len(A)), scan=None))
            continue
        scan = stability_scan(R, A, caps=caps, maxdim=3)
        rec = dict(tau=float(tc), nA=int(len(A)),
                   scan={int(c): (list(v[0]) if v[0] != "too_big" else "too_big")
                         for c, v in scan.items()})
        if check_pm:
            layer = layer_index(R, A, nmax=5)
            pm = {}
            for c in caps:
                mx = mrs_nerve_fast(R, A, layer, 4, shadow_cap=c)
                if max((len(s) for s in mx), default=0) <= 13:
                    pm[int(c)] = bool(is_pseudomanifold(mx))
            rec["pseudomanifold"] = pm
        out.append(rec)

    seq = []
    for s in out:
        if s["scan"] is None:
            seq.append(None); continue
        vals = [tuple(v) for v in s["scan"].values() if v != "too_big"]
        seq.append(max(set(vals), key=vals.count) if vals else None)
    chi, events = morse_count([v for v in seq if v is not None])

    return dict(N=N, seed=int(seed), a0=a0, t_build=t_build,
                t_total=time.time() - t0, slices=out,
                plateau=[list(v) if v else None for v in seq],
                morse_chi=int(chi),
                events=[[int(a), int(b)] for a, b in events])


def scaling_study(seeds=3, caps=(6, 8, 10, 12), Ns=(3000, 10000, 30000)):
    print("=" * 76)
    print(" SCALING STUDY -- run this BEFORE requesting a large allocation")
    print("=" * 76)
    print(f"\n{'N':>8} {'|A|':>6} {'plateau width':>15} {'modal beta':>16} {'time':>8}")
    for N in Ns:
        widths, modal, tt = [], [], []
        for s in range(seeds):
            t0 = time.time()
            tau, x, R, meta = gs_wormhole(N, seed=s)
            A = antichain_from_band(R, tau, 0.0, 0.10)
            if len(A) < 10:
                continue
            scan = stability_scan(R, A, caps=caps, maxdim=3)
            good = [tuple(v[0]) for v in scan.values() if v[0] != "too_big"]
            widths.append(len(good))
            if good:
                modal.append(max(set(good), key=good.count))
            tt.append(time.time() - t0)
            nA = len(A)
        if not tt:
            print(f"{N:>8} {'--':>6}")
            continue
        mm = max(set(modal), key=modal.count) if modal else None
        print(f"{N:>8} {nA:>6} {np.mean(widths):>15.1f} {str(mm):>16} "
              f"{np.mean(tt):>7.1f}s")
    print("\n Interpretation: if the plateau width DECREASES with N, the")
    print(" observable is not converging and a larger run will not help.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--seed0", type=int, default=0)
    ap.add_argument("--nslices", type=int, default=9)
    ap.add_argument("--band", type=float, default=0.10)
    ap.add_argument("--a0", type=float, default=1.0)
    ap.add_argument("--caps", type=int, nargs="+", default=[6, 8, 10, 12])
    ap.add_argument("--pseudomanifold", action="store_true")
    ap.add_argument("--scaling", action="store_true")
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--mpi", action="store_true")
    args = ap.parse_args()

    if args.scaling:
        scaling_study()
        return

    rank, size = args.task, args.ntasks
    if args.mpi:
        from mpi4py import MPI                      # noqa: F401
        comm = MPI.COMM_WORLD
        rank, size = comm.Get_rank(), comm.Get_size()

    seeds = [args.seed0 + s for s in range(args.seeds)][rank::size]
    recs = []
    for s in seeds:
        r = one_seed(args.N, s, tuple(args.caps), args.nslices, args.band,
                     a0=args.a0, check_pm=args.pseudomanifold)
        recs.append(r)
        print(f"[rank {rank}] seed {s}: chi={r['morse_chi']:+d} "
              f"({r['t_total']:.1f}s)", flush=True)

    if args.out:
        os.makedirs(args.out, exist_ok=True)
        json.dump(recs, open(os.path.join(args.out, f"task{rank:04d}.json"), "w"))
    else:
        from collections import Counter
        print("\nMorse chi:", dict(Counter(r["morse_chi"] for r in recs)))
        print("target: -2 for a wormhole (index 1 + index 3)")


if __name__ == "__main__":
    main()
