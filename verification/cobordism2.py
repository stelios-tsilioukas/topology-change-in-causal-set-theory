import numpy as np, gudhi
from homology import betti_and_chi
from collections import Counter

# ---------- (A) EXACT bulk chi: triangulated pair of pants = rectangle with 2 holes ----------
def triangulated_pants(nx=13, ny=6, holes=((3,1,2,2),(8,1,2,2))):
    def vid(i,j): return j*nx+i
    def in_hole(ci,cj):
        return any(i0<=ci<i0+w and j0<=cj<j0+h for (i0,j0,w,h) in holes)
    tris=[]
    for i in range(nx-1):
        for j in range(ny-1):
            if in_hole(i,j): continue
            a,b,c,d=vid(i,j),vid(i+1,j),vid(i,j+1),vid(i+1,j+1)
            tris.append(tuple(sorted((a,b,c)))); tris.append(tuple(sorted((b,d,c))))
    return tris

print("=== (A) EXACT bulk Euler characteristic of pair-of-pants cobordism W^2 ===")
tris=triangulated_pants()
b,chi,f=betti_and_chi(tris)
print(f"  triangulated pants: betti={[b[k] for k in sorted(b)]}  chi={chi}   (V,E,F)=({f.get(0)},{f.get(1)},{f.get(2)})")
print(f"  => bulk chi = {chi}   [pair of pants: b=[1,2,0], chi=-1]")
print(f"  boundary = 3 circles (1 'incoming' waist + 2 'outgoing' legs); each boundary circle has chi=0\n")

# ---------- (B) Causal-set spatial slicing: 1 circle -> 2 circles ----------
def sample_tube(A,B,rho,m,rng):
    A=np.array(A,float);B=np.array(B,float);ax=B-A;L=np.linalg.norm(ax);ax/=L
    tmp=np.array([1,0,0.]) if abs(ax[0])<0.9 else np.array([0,1,0.])
    u=np.cross(ax,tmp);u/=np.linalg.norm(u);v=np.cross(ax,u)
    s=rng.uniform(0,L,m);phi=rng.uniform(0,2*np.pi,m)
    return A[None,:]+s[:,None]*ax[None,:]+rho*(np.cos(phi)[:,None]*u[None,:]+np.sin(phi)[:,None]*v[None,:])

def make_pants(rng,rho=0.28,m=1500):
    C=np.array([0,0,0.]);Pw=np.array([0,0,1.3]);Pl=np.array([-0.85,0,-1.15]);Pr=np.array([0.85,0,-1.15])
    return np.vstack([sample_tube(C,Pw,rho,m,rng),sample_tube(C,Pl,rho,m,rng),sample_tube(C,Pr,rho,m,rng)])

def slice_topology(xy):
    """widest connected... here count components (b0) and loops (b1) via alpha persistence."""
    st=gudhi.AlphaComplex(points=xy).create_simplex_tree(); st.compute_persistence()
    def cnt(d):
        iv=st.persistence_intervals_in_dimension(d)
        fin=sorted([de-bi for (bi,de) in iv if not np.isinf(de)],reverse=True)
        ninf=sum(1 for (bi,de) in iv if np.isinf(de))
        # signal bars: lifetime > 6x median-ish noise; use gap detection
        sig=0
        if fin:
            thr=0.15*fin[0]
            sig=sum(1 for x in fin if x>thr)
        return ninf,sig
    b0inf,_=cnt(0); _,b1=cnt(1)
    return b0inf,b1

print("=== (B) Causal-set spatial slices (antichains via global time z) ===")
rng=np.random.default_rng(1); P=make_pants(rng,m=2500); z=P[:,2]
early=P[(z>0.5)&(z<1.0)][:,:2]
late =P[(z<-0.55)&(z>-1.05)][:,:2]
b0e,b1e=slice_topology(early); b0l,b1l=slice_topology(late)
print(f"  early antichain (high z, waist): #circles(b0)={b0e}, loops(b1)={b1e}")
print(f"  late  antichain (low z, legs)  : #circles(b0)={b0l}, loops(b1)={b1l}")
print(f"  spatial change: {b0e} circle -> {b0l} circles   (Delta components = +{b0l-b0e})")
print(f"  chi(spatial slice) = 0 for every slice (1-manifolds)\n")

print("=== (C) The Morse bridge ===")
print("  1 circle -> 2 circles is ONE index-1 saddle (a 'pair of pants' critical point).")
print("  Cobordism formula:  chi(W) = chi(Sigma_0) + sum_i (-1)^lambda_i = 0 + (-1)^1 = -1.")
print(f"  Independently measured bulk chi = {chi}.   AGREEMENT: {chi==-1}\n")

print("=== (D) The physical 3+1D bridge (analytic, Kunneth) ===")
for name,factors,val in [("S^1xS^3 (Euclid. wormhole)",[0,0],0),
                         ("S^2xS^2 (Nariai)",[2,2],4)]:
    prod=factors[0]*factors[1]
    dchi=prod-2
    print(f"  {name}: chi={prod}, connected-sum insert => dchi = chi-2 = {dchi:+d}")
print("  => matches the paper's wormhole dchi=-2 and Nariai dchi=+2.")
