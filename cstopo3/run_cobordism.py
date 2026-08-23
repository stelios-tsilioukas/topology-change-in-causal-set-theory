#!/usr/bin/env python3
"""
run_cobordism.py -- the 3+1D counterpart of the 2+1D trousers demonstration.

Detects the wormhole-sector topology change

    S^3  --(index-1 surgery)-->  S^1 x S^2 ,      Delta beta = (0,+1,+1,0)

from the causal order alone.  A single sprinkling contains the transition; the
Betti numbers of a spatial slice are extracted before and after t_split and
compared.

    python run_cobordism.py --N 20000 --seeds 8
    python run_cobordism.py --N 60000 --seeds 128 --nL 350 --cap 16 \
        --task $SLURM_ARRAY_TASK_ID --ntasks 32 --out results/
    python run_cobordism.py --merge results/

Note on scope: the before- and after-slice witness windows lie entirely on their
own side of t_split, so the observable never uses a cross-split relation.  What
is measured is the pair of slice topologies, which is exactly Delta beta.
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cstopo3.geom3 import S3Warped, Surgery, Spacetime
from cstopo3.extract import slice_topology


def one_seed(N, seed, nL, cap, lam, neck, T, tsplit, maxdim, verbose=False):
    sb = S3Warped(1.0, neck=neck)
    sa = Surgery(1.0, lam, neck=neck)
    t0 = time.time()
    st = Spacetime(sb, N, T=T, seed=seed, slice_after=sa, t_split=tsplit)
    # place each antichain so that its witness window stays on its own side
    t_before = 0.25 * tsplit
    t_after = tsplit + 0.25 * (T - tsplit)
    bb, ib = slice_topology(st, t_before, nL, cap, seed=seed,
                            volume=sb.volume, maxdim=maxdim, verbose=verbose)
    ba, ia = slice_topology(st, t_after, nL, cap, seed=seed,
                            volume=sa.volume, maxdim=maxdim, verbose=verbose)
    want_b = (1, 0) if maxdim == 1 else (1, 0, 0)
    want_a = (1, 1) if maxdim == 1 else (1, 1, 1)
    ok = (bb == want_b and ba == want_a)
    return dict(N=N, seed=seed, nL=nL, cap=cap, lam=lam, neck=neck,
                beta_before=list(bb) if bb else None,
                beta_after=list(ba) if ba else None,
                delta=[y - x for x, y in zip(bb, ba)] if bb and ba else None,
                detected=bool(ok), info_before=ib, info_after=ia,
                wall=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--seed0", type=int, default=1)
    ap.add_argument("--nL", type=int, default=250, help="landmark count")
    ap.add_argument("--cap", type=int, default=14, help="shadow-cardinality cap")
    ap.add_argument("--lam", type=float, default=0.40, help="throat radius / R")
    ap.add_argument("--neck", type=float, default=1.5,
                    help="tube length at each end; the handle 1-cycle has "
                         "length 2*neck+pi and must exceed the cover radius")
    ap.add_argument("--T", type=float, default=4.5)
    ap.add_argument("--tsplit", type=float, default=2.0)
    ap.add_argument("--maxdim", type=int, default=1, choices=(1, 2),
                    help="1 -> (b0,b1) in seconds; 2 -> (b0,b1,b2) in minutes")
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--procs", type=int, default=1, help="local worker processes")
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--merge", type=str, default=None)
    a = ap.parse_args()

    if a.merge:
        rows = []
        for fn in sorted(os.listdir(a.merge)):
            if fn.endswith(".json"):
                rows += json.load(open(os.path.join(a.merge, fn)))["rows"]
        det = sum(r["detected"] for r in rows)
        print(f"merged {len(rows)} runs from {a.merge}")
        print(f"Delta beta detected: {det}/{len(rows)}  "
              f"({100.0*det/max(len(rows),1):.1f}%)")
        from collections import Counter
        print("beta before:", Counter(tuple(r["beta_before"] or ()) for r in rows))
        print("beta after :", Counter(tuple(r["beta_after"] or ()) for r in rows))
        print(f"mean wall per seed: {np.mean([r['wall'] for r in rows]):.1f}s")
        return

    seeds = [a.seed0 + k for k in range(a.seeds)][a.task::a.ntasks]
    print(f"# cobordism S^3 -> S^1xS^2   N={a.N} |L|={a.nL} cap={a.cap} "
          f"lam={a.lam} neck={a.neck} maxdim={a.maxdim}")
    print(f"# handle 1-cycle length = {2*a.neck+np.pi:.2f}, "
          f"{len(seeds)} seeds on task {a.task}/{a.ntasks}")

    args = [(a.N, s, a.nL, a.cap, a.lam, a.neck, a.T, a.tsplit, a.maxdim)
            for s in seeds]
    if a.procs > 1:
        import multiprocessing as mp
        with mp.Pool(a.procs) as pool:
            rows = pool.starmap(one_seed, args)
    else:
        rows = [one_seed(*x, verbose=True) for x in args]

    for r in rows:
        print(f"  seed {r['seed']:4d}  before={r['beta_before']} "
              f"after={r['beta_after']}  delta={r['delta']}  "
              f"{'DETECTED' if r['detected'] else 'no'}  {r['wall']:.0f}s")
    det = sum(r["detected"] for r in rows)
    print(f"\n  detection rate: {det}/{len(rows)}")

    if a.out:
        os.makedirs(a.out, exist_ok=True)
        fn = os.path.join(a.out, f"cob_task{a.task:04d}.json")
        json.dump(dict(config=vars(a), rows=rows), open(fn, "w"), indent=1)
        print(f"  wrote {fn}")


if __name__ == "__main__":
    main()
