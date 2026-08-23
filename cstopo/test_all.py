"""Regression tests. Run: python test_all.py"""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cstopo.causet import sprinkle_box, build_causal_flat
from cstopo.topology import betti_gf2, euler_char, antichain_from_band, stability_scan
from cstopo.geometry import trousers_2p1, cylinder_2p1

fails = 0
def check(name, got, exp):
    global fails
    ok = got == exp
    fails += (not ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}: {got}" + ("" if ok else f" != {exp}"))

print("== causal matrix vs brute force ==")
t, x = sprinkle_box(300, 2.0, 2.0, 4, seed=1)
R = build_causal_flat(t, x)
dt = t[None,:]-t[:,None]; dx2 = ((x[None,:,:]-x[:,None,:])**2).sum(-1)
B = (dt>0)&(dt**2>dx2); np.fill_diagonal(B,False)
check("dense match", bool(np.array_equal(R.to_dense(), B)), True)

print("== homology engine ==")
check("circle", betti_gf2([(0,1),(1,2),(0,2)])[:2], (1,1))
check("S^2", betti_gf2([(a,b,c) for a in (0,1) for b in (2,3) for c in (4,5)])[:3], (1,0,1))
def vid(i,j): return 3*(i%3)+(j%3)
tor = list({tuple(sorted(u)) for i in range(3) for j in range(3) for u in
     [(vid(i,j),vid(i+1,j),vid(i,j+1)),(vid(i+1,j),vid(i+1,j+1),vid(i,j+1))]})
check("T^2", betti_gf2(tor)[:3], (1,2,1))
check("chi(T^2)", euler_char(tor), 0)

print("== order-only slice topology ==")
t,x,R,m = cylinder_2p1(2000, seed=0)
A = antichain_from_band(R, t, 1.0, 0.10)
check("cylinder", tuple(stability_scan(R,A,caps=(10,),maxdim=2)[10][0][:2]), (1,1))

t,x,R,m = trousers_2p1(6000, seed=0)
Aw = antichain_from_band(R, t, 0.5, 0.08)
Al = antichain_from_band(R, t, 1.7, 0.08)
check("trousers waist", tuple(stability_scan(R,Aw,caps=(10,),maxdim=2)[10][0][:2]), (1,1))
check("trousers legs ", tuple(stability_scan(R,Al,caps=(10,),maxdim=2)[10][0][:2]), (2,2))

print(f"\n{'ALL TESTS PASSED' if fails==0 else str(fails)+' FAILURES'}")
sys.exit(1 if fails else 0)
