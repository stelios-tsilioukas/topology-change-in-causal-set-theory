#!/usr/bin/env python3
"""
run_static.py -- static-slice validation and the resolution plateau.

Sprinkles ultrastatic cylinders Sigma x [0,T] with Sigma = S^3 and S^1 x S^2,
both with EXACT closed-form geodesic distance, and recovers beta from the causal
order alone.  This is the make-or-break test: if the two manifolds the wormhole
surgery interpolates cannot be told apart on a static slab, the cobordism is
hopeless.  It also locates the plateau in the shadow-cardinality resolution.

    python run_static.py --N 20000 --seeds 8 --caps 8 10 12 14 16
    python run_static.py --N 40000 --seeds 8 --nL 150 --caps 12 --maxdim 2
    python run_static.py --scaling
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cstopo3.geom3 import S3, S1xS2, S3Warped, Surgery, Spacetime
from cstopo3.extract import slice_topology

GEOMS = {
    "S3":       lambda: S3(1.0),
    "S1xS2":    lambda: S1xS2(1.0, 1.0),
    "S3warped": lambda: S3Warped(1.0, neck=1.5),
    "surgery":  lambda: Surgery(1.0, 0.40, neck=1.5),
}


def one(gname, N, seed, nL, cap, maxdim, T=3.0):
    sl = GEOMS[gname]()
    st = Spacetime(sl, N, T=T, seed=seed)
    b, info = slice_topology(st, 0.4 * T, nL, cap, seed=seed,
                             volume=sl.volume, maxdim=maxdim)
    want = sl.betti[:maxdim + 1]
    return dict(geom=gname, N=N, seed=seed, nL=nL, cap=cap,
                betti=list(b) if b else None, want=list(want),
                ok=bool(b == want), info=info)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geoms", nargs="+", default=["S3", "S1xS2"],
                    choices=list(GEOMS))
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--seed0", type=int, default=1)
    ap.add_argument("--nL", type=int, default=250)
    ap.add_argument("--caps", type=int, nargs="+", default=[8, 10, 12, 14, 16])
    ap.add_argument("--maxdim", type=int, default=1, choices=(1, 2))
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--scaling", action="store_true",
                    help="cost and |A| versus N; run this first on new hardware")
    a = ap.parse_args()

    if a.scaling:
        print("# scaling: |A|, wall time per slice, cap=14, |L|=250, maxdim=1")
        print(f"{'N':>8} {'geom':>8} {'ell':>8} {'|A|':>7} {'|L|':>5} "
              f"{'facets':>7} {'wall':>7}  beta")
        for N in (5000, 10000, 20000, 40000, 80000):
            for g in ("S3", "S1xS2"):
                nL = min(250, int(0.15 * N ** 0.75))
                t0 = time.time()
                try:
                    r = one(g, N, 1, nL, 14, 1)
                except Exception as e:
                    print(f"{N:>8} {g:>8}  FAILED: {e}"); continue
                i = r["info"]
                print(f"{N:>8} {g:>8} {i.get('ell',0):8.4f} {i.get('nA',0):>7} "
                      f"{i.get('nL',0):>5} {i.get('n_facets',0):>7} "
                      f"{time.time()-t0:6.1f}s  {r['betti']} "
                      f"{'OK' if r['ok'] else 'MISMATCH'}")
        return

    seeds = [a.seed0 + k for k in range(a.seeds)][a.task::a.ntasks]
    jobs = [(g, a.N, s, a.nL, c, a.maxdim)
            for g in a.geoms for c in a.caps for s in seeds]
    print(f"# static  N={a.N} |L|={a.nL} caps={a.caps} maxdim={a.maxdim} "
          f"seeds={len(seeds)}  ({len(jobs)} runs)")

    if a.procs > 1:
        import multiprocessing as mp
        with mp.Pool(a.procs) as pool:
            rows = pool.starmap(one, jobs)
    else:
        rows = [one(*j) for j in jobs]

    print(f"\n{'geom':>9} {'cap':>4}  detection   modal beta")
    for g in a.geoms:
        for c in a.caps:
            sub = [r for r in rows if r["geom"] == g and r["cap"] == c]
            if not sub:
                continue
            from collections import Counter
            modal = Counter(tuple(r["betti"] or ()) for r in sub).most_common(1)[0]
            hits = sum(r["ok"] for r in sub)
            print(f"{g:>9} {c:>4}   {hits}/{len(sub)}       {modal[0]} "
                  f"(x{modal[1]})   want {sub[0]['want']}")

    if a.out:
        os.makedirs(a.out, exist_ok=True)
        fn = os.path.join(a.out, f"static_task{a.task:04d}.json")
        json.dump(dict(config=vars(a), rows=rows), open(fn, "w"), indent=1)
        print(f"  wrote {fn}")


if __name__ == "__main__":
    main()
