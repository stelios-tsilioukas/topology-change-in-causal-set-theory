import numpy as np, sys, time; sys.path.insert(0,'.')
from cstopo.geometry import trousers_2p1, cylinder_2p1
from cstopo.topology import antichain_from_band, stability_scan
CAPS=(4,6,8,10,12,14)
print("=== cylinder control, N=8000 ===")
t,x,R,m=cylinder_2p1(8000,seed=0); A=antichain_from_band(R,t,1.0,0.08)
sc=stability_scan(R,A,caps=CAPS,maxdim=2)
print(f" |A|={len(A)}: "+" ".join(f"{c}:{v[0][:2] if v[0]!='too_big' else 'big'}" for c,v in sc.items()))
print("\n=== trousers N=8000, full cap sweep ===")
t,x,R,m=trousers_2p1(8000,seed=0)
for lab,tc in [("waist",0.5),("legs",1.7)]:
    A=antichain_from_band(R,t,tc,0.08); sc=stability_scan(R,A,caps=CAPS,maxdim=2)
    print(f" {lab:>6} |A|={len(A):>4}: "+" ".join(f"{c}:{v[0][:2] if v[0]!='too_big' else 'big'}" for c,v in sc.items()))
print("\n=== timing & antichain size vs N ===")
print(f"{'N':>7} {'|A|':>6} {'t_build':>9} {'t_topo':>9}")
for N in (3000,6000,12000,24000):
    t0=time.time(); t,x,R,m=trousers_2p1(N,seed=2); tb=time.time()-t0
    t1=time.time(); A=antichain_from_band(R,t,1.7,0.08)
    stability_scan(R,A,caps=(10,),maxdim=2); tt=time.time()-t1
    print(f"{N:>7} {len(A):>6} {tb:>8.1f}s {tt:>8.1f}s")
