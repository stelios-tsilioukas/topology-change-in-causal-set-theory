#!/usr/bin/env python3
"""
run_b2.py -- convergence of b_2, the expensive observable.  THE CLUSTER JOB.

b_0 and b_1 are cheap: the boundary of a triangle, written in the fundamental
cycle basis of a spanning forest, has at most three nonzero coordinates, so the
rank computation is trivial (see cstopo3/fast01.py).  b_2 has no such shortcut
-- it needs the 3-simplices, C(cap,4) of them per facet -- and it converges only
at larger cap, because filling a 3-dimensional void requires the complex to be
solid in dimension 3.

Measured on a single sandbox core, S^3, N = 20000, |L| = 250:

    cap    b2     raw 3-simplices    wall
     10    42          68,270          5 s
     12    26         171,395         19 s
     14     6         383,321         67 s
     16     3         710,921        198 s

The trend is clean and heads to the correct b2 = 0; cap ~ 18-20 should reach it,
at an extrapolated 550-1400 s per (geometry, cap, seed).  That is what this
script is for.  b_1 is unaffected and correct throughout.

    # one array task per (geometry, cap, seed) triple
    python run_b2.py --caps 14 16 18 20 --seeds 8 \
        --task $SLURM_ARRAY_TASK_ID --ntasks 32 --out results_b2/
    python run_b2.py --merge results_b2/
"""
from __future__ import annotations
import argparse, json, os, sys, time
import numpy as np
from math import comb

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cstopo3.geom3 import S3, S1xS2, S3Warped, Surgery, Spacetime
from cstopo3.extract import antichain, landmarks, landmark_shadows
from cstopo3.z2 import betti, faces_by_dim

GEOMS = {
    "S3":       lambda: S3(1.0),
    "S1xS2":    lambda: S1xS2(1.0, 1.0),
    "S3warped": lambda: S3Warped(1.0, neck=1.5),
    "surgery":  lambda: Surgery(1.0, 0.40, neck=1.5),
}


def one(gname, N, seed, nL, cap, T=3.0, guard=40_000_000):
    sl = GEOMS[gname]()
    t0 = time.time()
    st = Spacetime(sl, N, T=T, seed=seed)
    A = antichain(st, t0=0.4 * T)
    if A.size < nL:
        return dict(geom=gname, N=N, seed=seed, nL=nL, cap=cap,
                    betti=None, note=f"|A|={A.size}<nL")
    L = landmarks(A, nL, seed=seed)
    tA = float(st.t[A].mean())
    facets, nW, h, cov = landmark_shadows(st, L, cap, tA, sl.volume)
    raw3 = sum(comb(len(f), 4) for f in facets if len(f) >= 4)
    if raw3 > guard:
        return dict(geom=gname, N=N, seed=seed, nL=nL, cap=cap, betti=None,
                    note=f"raw3={raw3} exceeds guard {guard}", raw3=raw3)
    b = betti(facets, maxdim=2)
    return dict(geom=gname, N=N, seed=seed, nL=nL, cap=cap,
                betti=list(b), want=list(sl.betti[:3]),
                ok=bool(tuple(b) == sl.betti[:3]),
                nA=int(A.size), nL_used=int(L.size), n_facets=len(facets),
                coverage=cov, height=float(h), raw3=int(raw3),
                wall=time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--geoms", nargs="+", default=["S3", "S1xS2"], choices=list(GEOMS))
    ap.add_argument("--N", type=int, default=20000)
    ap.add_argument("--nL", type=int, default=250)
    ap.add_argument("--caps", type=int, nargs="+", default=[14, 16, 18, 20])
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--seed0", type=int, default=1)
    ap.add_argument("--procs", type=int, default=1)
    ap.add_argument("--task", type=int, default=0)
    ap.add_argument("--ntasks", type=int, default=1)
    ap.add_argument("--guard", type=int, default=40_000_000,
                    help="skip if raw 3-simplex count exceeds this (memory)")
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--merge", type=str, default=None)
    a = ap.parse_args()

    if a.merge:
        rows = []
        for fn in sorted(os.listdir(a.merge)):
            if fn.startswith("b2_") and fn.endswith(".json"):
                rows += json.load(open(os.path.join(a.merge, fn)))["rows"]
        print(f"merged {len(rows)} runs\n")
        print(f"{'geom':>9} {'cap':>4} {'n':>3} {'b2 values':>28} {'ok':>6} {'mean wall':>10}")
        for g in sorted({r["geom"] for r in rows}):
            for c in sorted({r["cap"] for r in rows if r["geom"] == g}):
                sub = [r for r in rows if r["geom"] == g and r["cap"] == c
                       and r.get("betti")]
                if not sub:
                    continue
                b2 = [r["betti"][2] for r in sub]
                ok = sum(r.get("ok", False) for r in sub)
                print(f"{g:>9} {c:>4} {len(sub):>3} {str(sorted(b2)):>28} "
                      f"{ok}/{len(sub):<4} {np.mean([r['wall'] for r in sub]):9.0f}s")
        return

    jobs = [(g, a.N, a.seed0 + k, a.nL, c)
            for g in a.geoms for c in a.caps for k in range(a.seeds)]
    jobs = jobs[a.task::a.ntasks]
    print(f"# b2 convergence: {len(jobs)} runs on task {a.task}/{a.ntasks}")

    if a.procs > 1:
        import multiprocessing as mp
        with mp.Pool(a.procs) as pool:
            rows = pool.starmap(one, jobs)
    else:
        rows = []
        for j in jobs:
            r = one(*j)
            rows.append(r)
            print(f"  {r['geom']:>9} cap={r['cap']:>3} seed={r['seed']:>3} "
                  f"beta={r.get('betti')} raw3={r.get('raw3','-')} "
                  f"{r.get('wall',0):.0f}s {r.get('note','')}")

    if a.out:
        os.makedirs(a.out, exist_ok=True)
        tag = f"N{a.N}_L{a.nL}_" + "-".join(str(c) for c in a.caps)
        fn = os.path.join(a.out, f"b2_{tag}_task{a.task:04d}.json")
        json.dump(dict(config=vars(a), rows=rows), open(fn, "w"), indent=1)
        print(f"  wrote {fn}")


if __name__ == "__main__":
    main()
