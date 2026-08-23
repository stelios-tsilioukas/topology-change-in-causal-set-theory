import sys, os, json, time; sys.path.insert(0,'.')
import numpy as np
import multiprocessing as mp
from cstopo3.geom3 import S3, S1xS2, Spacetime
from cstopo3.extract import antichain, landmarks, landmark_shadows
from cstopo3.fast01 import betti01
from cstopo3.z2 import betti

GE = {"S3": lambda: S3(1.0), "S1xS2": lambda: S1xS2(1.0,1.0)}

def plateau_job(args):
    g, seed, caps, N, nL = args
    sl = GE[g]()
    st = Spacetime(sl, N, T=3.0, seed=seed)
    A = antichain(st, t0=1.2); L = landmarks(A, nL, seed=seed)
    tA = float(st.t[A].mean()); out={}
    for c in caps:
        f,nW,h,cov = landmark_shadows(st, L, c, tA, sl.volume)
        if not f: out[c]=None; continue
        (b0,b1),_ = betti01(f, L.size)
        out[c]=[b0,b1,cov,len(f)]
    return dict(geom=g, seed=seed, N=N, nL=nL, res={str(k):v for k,v in out.items()})

def b2_job(args):
    g, seed, c, N, nL = args
    sl = GE[g]()
    t0=time.time()
    st = Spacetime(sl, N, T=3.0, seed=seed)
    A = antichain(st, t0=1.2); L = landmarks(A, nL, seed=seed)
    f,nW,h,cov = landmark_shadows(st, L, c, float(st.t[A].mean()), sl.volume)
    b = betti(f, maxdim=2)
    return dict(geom=g, seed=seed, cap=c, betti=list(b), wall=time.time()-t0)

if __name__ == "__main__":
    caps=[6,8,10,12,14,16,18]
    jobs=[(g,s,caps,20000,250) for g in GE for s in range(1,9)]
    with mp.Pool(2) as p: plate = p.map(plateau_job, jobs)
    json.dump(plate, open("fig_plateau.json","w"), indent=1)
    print("plateau done", flush=True)

    j2=[(g,s,c,20000,250) for g in GE for c in (10,12,14) for s in (1,2,3)]
    j2+=[(g,1,16,20000,250) for g in GE]
    with mp.Pool(2) as p: b2 = p.map(b2_job, j2)
    json.dump(b2, open("fig_b2.json","w"), indent=1)
    print("b2 done", flush=True)
