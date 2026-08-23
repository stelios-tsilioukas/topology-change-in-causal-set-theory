import numpy as np, gudhi
from collections import Counter

def sphere_pts(n,seed):
    rng=np.random.default_rng(seed); v=rng.normal(size=(n,3)); v/=np.linalg.norm(v,axis=1,keepdims=True); return v
def torus_pts(n,seed,Rmaj=2.0,rmin=0.8):
    rng=np.random.default_rng(seed)
    th=rng.uniform(0,2*np.pi,n); ph=rng.uniform(0,2*np.pi,n)
    return np.column_stack([(Rmaj+rmin*np.cos(th))*np.cos(ph),(Rmaj+rmin*np.cos(th))*np.sin(ph),rmin*np.sin(th)])

def intervals(points):
    st=gudhi.AlphaComplex(points=points).create_simplex_tree()
    st.compute_persistence()
    return {d:st.persistence_intervals_in_dimension(d) for d in (0,1,2)}

def betti_at(iv,f):
    b={}
    for d in (0,1,2):
        b[d]=sum(1 for (bi,de) in iv[d] if bi<=f and (np.isinf(de) or f<de))
    return (b[0],b[1],b[2])

def widest_plateau(points,ngrid=400):
    iv=intervals(points)
    fmax=max((de for d in (0,1,2) for (bi,de) in iv[d] if not np.isinf(de)),default=1.0)
    grid=np.linspace(0.0,fmax,ngrid)
    vecs=[betti_at(iv,f) for f in grid]
    # widest contiguous run of an identical betti vector, restricted to connected (b0==1)
    best=None; best_len=-1; i=0
    while i<len(vecs):
        j=i
        while j+1<len(vecs) and vecs[j+1]==vecs[i]: j+=1
        length=grid[j]-grid[i]
        if vecs[i][0]==1 and length>best_len:
            best_len=length; best=vecs[i]
        i=j+1
    return best

def study(name,ptfun,n,seeds=12):
    rows=[widest_plateau(ptfun(n,s)) for s in range(seeds)]
    modal=Counter(rows).most_common(3)
    print(f"{name} (n={n}, {seeds} samples): widest connected Betti plateau")
    for bt,c in modal:
        print(f"    betti={list(bt)}  chi={bt[0]-bt[1]+bt[2]:+d}   [{c}/{seeds}]")
    top=modal[0][0]
    print(f"    -> modal chi = {top[0]-top[1]+top[2]:+d}\n")
    return top

if __name__=="__main__":
    print("=== Widest-plateau Betti (rigorous MRS-style stability window) ===\n")
    bS=study("S^2 sphere", sphere_pts, 500)
    bT=study("T^2 torus",  torus_pts,  700)
    chiS=bS[0]-bS[1]+bS[2]; chiT=bT[0]-bT[1]+bT[2]
    print(f">>> S^2: betti={list(bS)}, chi={chiS:+d}")
    print(f">>> T^2: betti={list(bT)}, chi={chiT:+d}")
    print(f">>> delta_chi (handle birth) = {chiT-chiS:+d}   [target: -2]")
