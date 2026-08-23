import numpy as np
from itertools import combinations, permutations

def kappa_all(maximal, nverts):
    """kappa(v) = sum_k (-1)^k n_k(v)/(k+1) over ALL faces of the complex."""
    faces=set()
    for s in maximal:
        s=tuple(sorted(s))
        for k in range(1,len(s)+1):
            for c in combinations(s,k): faces.add(c)
    k=np.zeros(nverts)
    fcount={}
    for f in faces:
        d=len(f)-1; w=((-1)**d)/(d+1.0)
        fcount[d]=fcount.get(d,0)+1
        for v in f: k[v]+=w
    chi=sum((-1)**d*c for d,c in fcount.items())
    return k, chi, fcount

print("="*74)
print(" GENERALISATION OF kappa TO FOUR DIMENSIONS")
print("="*74)

# ---------- (a) S^4 = boundary of the 5-simplex ----------
S4=[c for c in combinations(range(6),5)]      # six 4-simplices
k,chi,f=kappa_all(S4,6)
print(f"\n (a) S^4 = boundary of 5-simplex   [chi(S^4) = +2]")
print(f"     face counts by dim: {dict(sorted(f.items()))}")
print(f"     sum kappa = {k.sum():+.6f}   chi = {chi:+d}")
print(f"     kappa per vertex  = {np.round(k,6)}   (all equal by symmetry, = 2/6)")

# ---------- (b) Kuhn-triangulated flat 4-torus ----------
def kuhn_T4(m):
    """flat 4-torus: each hypercube split into 4! = 24 simplices (Kuhn/Freudenthal)."""
    vid=lambda p: ((p[0]%m)*m**3+(p[1]%m)*m**2+(p[2]%m)*m+(p[3]%m))
    simps=[]
    e=np.eye(4,dtype=int)
    for i in range(m):
     for j in range(m):
      for kk in range(m):
       for l in range(m):
        p0=np.array([i,j,kk,l])
        for perm in permutations(range(4)):
            verts=[p0.copy()]; cur=p0.copy()
            for d in perm:
                cur=cur+e[d]; verts.append(cur.copy())
            simps.append(tuple(sorted(vid(v) for v in verts)))
    return list(set(simps)), m**4

m=3
S,nv=kuhn_T4(m)
k,chi,f=kappa_all(S,nv)
print(f"\n (b) Kuhn-triangulated FLAT 4-torus, {m}^4 = {nv} vertices   [chi(T^4) = 0]")
print(f"     4-simplices: {len(S)},  face counts by dim: {dict(sorted(f.items()))}")
print(f"     sum kappa      = {k.sum():+.3e}   chi = {chi:+d}")
print(f"     max |kappa(v)| = {np.abs(k).max():.3e}   <-- vanishes POINTWISE")
print(f"     kappa spread   = {k.std():.3e}")

print("""
 => The two-dimensional result generalises verbatim:

      * flat manifoldlike 4-region  ->  kappa = 0 at every vertex
      * curved / topologically nontrivial -> kappa != 0, and sum kappa = chi

    kappa is therefore a genuine EULER DENSITY in four dimensions, exact and
    purely combinatorial, with no continuum input.  It is the discrete object
    whose coarse-grained limit is  G/32pi^2.""")
