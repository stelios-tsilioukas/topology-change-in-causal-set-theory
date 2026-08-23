import numpy as np
from kappa2 import kappa_from_simplices

def torus_grid(M,N):
    """standard triangulated flat torus: every vertex degree 6."""
    vid=lambda i,j: (i%M)*N+(j%N)
    tris=set()
    for i in range(M):
        for j in range(N):
            a,b,c,d=vid(i,j),vid(i+1,j),vid(i,j+1),vid(i+1,j+1)
            tris.add(tuple(sorted((a,b,c)))); tris.add(tuple(sorted((b,d,c))))
    return list(tris), M*N

def all_simplices(tris):
    S=set()
    for t in tris:
        S.add(t)
        for e in ((t[0],t[1]),(t[0],t[2]),(t[1],t[2])): S.add(e)
        for v in t: S.add((v,))
    return list(S)

def degrees(tris,nv):
    deg=np.zeros(nv,int)
    E=set()
    for t in tris:
        for e in ((t[0],t[1]),(t[0],t[2]),(t[1],t[2])): E.add(e)
    for a,b in E: deg[a]+=1; deg[b]+=1
    return deg

M,N=6,6
tris,nv=torus_grid(M,N)
k=kappa_from_simplices(all_simplices(tris),nv)
deg=degrees(tris,nv)
print("="*74)
print(" (A) FLAT TRIANGULATED TORUS  (6x6, every vertex degree 6)")
print("="*74)
print(f"   degrees: all = {set(deg.tolist())}")
print(f"   kappa:   max|kappa| = {np.abs(k).max():.2e}   sum = {k.sum():+.6f}   chi(T^2)=0")
print("   => the Euler density vanishes POINTWISE on a flat manifoldlike region.")

# ---- introduce ONE topological defect by an edge flip ----
def flip(tris, M, N, i=2, j=2):
    vid=lambda a,b: (a%M)*N+(b%N)
    a,b,c,d=vid(i,j),vid(i+1,j),vid(i,j+1),vid(i+1,j+1)
    t1=tuple(sorted((a,b,c))); t2=tuple(sorted((b,d,c)))
    out=[t for t in tris if t!=t1 and t!=t2]
    out.append(tuple(sorted((a,b,d)))); out.append(tuple(sorted((a,d,c))))  # flipped diagonal
    return out, (a,d,b,c)

tris2,(a,d,b,c)=flip(tris,M,N)
k2=kappa_from_simplices(all_simplices(tris2),nv)
deg2=degrees(tris2,nv)
print("\n"+"="*74)
print(" (B) THE SAME TORUS WITH ONE EDGE FLIP  (a single topological defect)")
print("="*74)
nz=np.nonzero(np.abs(k2)>1e-9)[0]
print(f"   vertices with kappa != 0 :  {nz.tolist()}")
print(f"   their degrees            :  {deg2[nz].tolist()}")
print(f"   their kappa              :  {np.round(k2[nz],6).tolist()}")
print(f"   everywhere else          :  kappa = 0  (max |kappa| off-defect = "
      f"{np.abs(np.delete(k2,nz)).max():.2e})")
print(f"   sum kappa = {k2.sum():+.6f}    (chi still 0: the defect is a dipole)")
print(f"\n   kappa(deg 5) = 1 - 5/6 = {1-5/6:+.6f}   kappa(deg 7) = 1 - 7/6 = {1-7/6:+.6f}")

print("\n"+"="*74)
print(" INTERPRETATION")
print("="*74)
print("""   kappa(v) = 1 - deg(v)/6  for a triangulated surface: EXACTLY the Regge
   angle deficit.  Hence

     * flat manifoldlike region        ->  kappa = 0 pointwise
     * topological defect              ->  kappa != 0, supported ON the defect
     * sum over all vertices           ->  chi   (exact, always)

   This is the atomic Euler measure  mu_chi = sum_i dchi_i delta(x-x_i)
   DERIVED rather than posited: the Euler density of a discrete space is
   automatically a sum of delta functions on its defects, and the smooth
   Gauss-Bonnet density is what a coarse-grained cloud of them looks like.""")
