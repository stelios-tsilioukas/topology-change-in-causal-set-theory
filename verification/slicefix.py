import numpy as np, gudhi
def sample_tube(A,B,rho,m,rng):
    A=np.array(A,float);B=np.array(B,float);ax=B-A;L=np.linalg.norm(ax);ax/=L
    tmp=np.array([1,0,0.]) if abs(ax[0])<0.9 else np.array([0,1,0.])
    u=np.cross(ax,tmp);u/=np.linalg.norm(u);v=np.cross(ax,u)
    s=rng.uniform(0,L,m);phi=rng.uniform(0,2*np.pi,m)
    return A[None,:]+s[:,None]*ax[None,:]+rho*(np.cos(phi)[:,None]*u[None,:]+np.sin(phi)[:,None]*v[None,:])
def make_pants(rng,rho=0.24,m=2500):
    C=np.array([0,0,0.]);Pw=np.array([0,0,1.3]);Pl=np.array([-0.95,0,-1.15]);Pr=np.array([0.95,0,-1.15])
    return np.vstack([sample_tube(C,Pw,rho,m,rng),sample_tube(C,Pl,rho,m,rng),sample_tube(C,Pr,rho,m,rng)])

def betti_at(iv,f,dims=(0,1)):
    return tuple(sum(1 for (bi,de) in iv[d] if bi<=f and (np.isinf(de) or f<de)) for d in dims)
def widest_plateau_2d(xy,ngrid=350):
    st=gudhi.AlphaComplex(points=xy).create_simplex_tree(); st.compute_persistence()
    iv={d:st.persistence_intervals_in_dimension(d) for d in (0,1)}
    fmax=max((de for d in (0,1) for (bi,de) in iv[d] if not np.isinf(de)),default=1.0)
    grid=np.linspace(0,fmax,ngrid); vecs=[betti_at(iv,f) for f in grid]
    best=None;bl=-1;i=0
    while i<len(vecs):
        j=i
        while j+1<len(vecs) and vecs[j+1]==vecs[i]: j+=1
        L=grid[j]-grid[i]
        # require at least one loop (a real slice), take widest such plateau
        if vecs[i][1]>=1 and L>bl: bl=L;best=vecs[i]
        i=j+1
    return best

rng=np.random.default_rng(1); P=make_pants(rng); z=P[:,2]
early=P[(z>0.5)&(z<1.05)][:,:2]
late =P[(z<-0.6)&(z>-1.1)][:,:2]
be=widest_plateau_2d(early); bl=widest_plateau_2d(late)
print("=== Causal-set spatial slices (widest-plateau Betti) ===")
print(f"  early antichain (waist): (b0,b1)={be}  -> {be[0]} circle")
print(f"  late  antichain (legs) : (b0,b1)={bl}  -> {bl[0]} circles")
print(f"  spatial topology change: {be[0]} -> {bl[0]} circles   (Delta components = +{bl[0]-be[0]})")
