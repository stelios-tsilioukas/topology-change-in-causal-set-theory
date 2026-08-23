import numpy as np, gudhi
from collections import Counter

def torus_pts(n,seed,R=2.0,r=0.8):
    rng=np.random.default_rng(seed)
    th=rng.uniform(0,2*np.pi,n); ph=rng.uniform(0,2*np.pi,n)
    return np.column_stack([(R+r*np.cos(th))*np.cos(ph),(R+r*np.cos(th))*np.sin(ph),r*np.sin(th)])

def widest_betti(P,ngrid=200):
    st=gudhi.AlphaComplex(points=P).create_simplex_tree(); st.compute_persistence()
    iv={d:st.persistence_intervals_in_dimension(d) for d in (0,1,2)}
    fmax=max((de for d in (0,1,2) for (b,de) in iv[d] if not np.isinf(de)),default=1.)
    grid=np.linspace(0,fmax,ngrid)
    vecs=[tuple(sum(1 for (b,de) in iv[d] if b<=f and (np.isinf(de) or f<de)) for d in (0,1,2)) for f in grid]
    best=None;bl=-1;i=0
    while i<len(vecs):
        j=i
        while j+1<len(vecs) and vecs[j+1]==vecs[i]: j+=1
        L=grid[j]-grid[i]
        if vecs[i][0]==1 and L>bl: bl=L;best=vecs[i]
        i=j+1
    return best

print("="*76)
print(" IS THE COARSE-GRAINING FLOW A POWER LAW OR A THRESHOLD?")
print("="*76)
print(" Sprinkle a torus (b_1 = 2), decimate at retention f, ask whether the")
print(" handle is still DETECTED.  Power law => gradual;  size effect => sharp.\n")
N0=1200
print(f"{'f':>7} {'N_kept':>8} {'P(detect b1>=1)':>18} {'P(b1=2)':>10}")
rows=[]
for f in (1.0,0.6,0.35,0.2,0.12,0.07,0.04,0.025,0.015):
    det=0; full=0; T=14
    for s in range(T):
        P=torus_pts(N0,s)
        rng=np.random.default_rng(1000+s)
        keep=rng.random(N0)<f
        Q=P[keep]
        if len(Q)<12: continue
        b=widest_betti(Q)
        if b is None: continue
        if b[1]>=1: det+=1
        if b[1]==2: full+=1
    rows.append((f,int(N0*f),det/T,full/T))
    print(f"{f:>7.3f} {int(N0*f):>8} {det/T:>18.2f} {full/T:>10.2f}")

print("""
 A power-law flow p_h -> f^(n_h-1) p_h would give a straight line of slope
 (n_h - 1) on a log-log plot of detection probability vs f.  Fit:""")
import numpy as np
xs=np.array([r[0] for r in rows]); ys=np.array([r[2] for r in rows])
m=(ys>0)&(ys<1)
if m.sum()>=3:
    sl,ic=np.polyfit(np.log(xs[m]),np.log(ys[m]),1)
    print(f"   slope over the transition region = {sl:.2f}  =>  effective n_h - 1 = {sl:.2f}")
    print(f"   i.e. effective n_h = {sl+1:.2f}")
else:
    print("   transition too sharp to fit a power law")
