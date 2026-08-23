import numpy as np
from itertools import combinations
from homology import build_complex, betti_and_chi

def kappa(maximal, nverts):
    """local combinatorial Gauss-Bonnet density
       kappa(v) = sum_k (-1)^k n_k(v)/(k+1),   n_k(v)=# k-simplices containing v"""
    by_dim = build_complex(maximal)
    k = np.zeros(nverts)
    for d, simps in by_dim.items():
        w = ((-1)**d)/(d+1.0)
        for s in simps:
            for v in s: k[v] += w
    return k

print("="*72)
print(" (1) IDENTITY CHECK:  chi = sum_v kappa(v)")
print("="*72)
# octahedron = S^2
oct_f=[(a,b,c) for a in (0,1) for b in (2,3) for c in (4,5)]
# 3x3 flat torus
def vid(i,j): return 3*(i%3)+(j%3)
tor_f=list({tuple(sorted(t)) for i in range(3) for j in range(3) for t in
     [(vid(i,j),vid(i+1,j),vid(i,j+1)),(vid(i+1,j),vid(i+1,j+1),vid(i,j+1))]})
tests=[("octahedron (S^2)",oct_f,6,2),("3x3 torus (T^2)",tor_f,9,0),
       ("filled triangle (disk)",[(0,1,2)],3,1)]
for name,f,nv,truechi in tests:
    k=kappa(f,nv); b,chi,_=betti_and_chi(f)
    print(f"  {name:>24}: sum kappa = {k.sum():+.6f}   chi = {chi:+d}   true = {truechi:+d}")
    print(f"  {'':>24}  kappa per vertex: {np.round(k,4)}")

print("\n"+"="*72)
print(" (2) IS kappa A LOCAL CURVATURE DENSITY?   (flat -> 0 pointwise)")
print("="*72)
print("  octahedron: all deg 4  -> kappa = 1 - 4/2 + 4/3 = 1/3 at every vertex (curved)")
print("  3x3 torus : all deg 6  -> kappa = 1 - 6/2 + 6/3 = 0   at every vertex (FLAT)")
print("  => kappa is exactly the Regge angle deficit: curvature concentrated at defects.\n")

# ---------- alpha-complex samples of S^2 and T^2 ----------
import gudhi
def alpha_maximal(points, alpha_max):
    st=gudhi.AlphaComplex(points=points).create_simplex_tree(max_alpha_square=alpha_max)
    simps=[tuple(s) for s,_ in st.get_simplices()]
    keep=[]
    S=set(simps)
    for s in simps:
        if not any(set(s)<set(t) for t in S): keep.append(s)
    return keep

def sphere_pts(n,seed):
    rng=np.random.default_rng(seed); v=rng.normal(size=(n,3))
    return v/np.linalg.norm(v,axis=1,keepdims=True)
def torus_pts(n,seed,R=2.0,r=0.8):
    rng=np.random.default_rng(seed)
    th=rng.uniform(0,2*np.pi,n); ph=rng.uniform(0,2*np.pi,n)
    return np.column_stack([(R+r*np.cos(th))*np.cos(ph),(R+r*np.cos(th))*np.sin(ph),r*np.sin(th)])

print("="*72)
print(" (3) kappa ON ALPHA COMPLEXES OF SAMPLED SURFACES")
print("="*72)
for name,fn,n,am,truechi in [("S^2",sphere_pts,220,0.12,2),("T^2",torus_pts,420,0.30,0)]:
    tot=[];sd=[]
    for s in range(4):
        P=fn(n,s); mx=alpha_maximal(P,am)
        k=kappa(mx,len(P)); tot.append(k.sum()); sd.append(k.std())
    print(f"  {name}:  sum kappa = {np.mean(tot):+.3f} +/- {np.std(tot):.3f}   (true chi = {truechi})")
    print(f"        pointwise spread std(kappa) = {np.mean(sd):.4f}")
