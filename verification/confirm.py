import numpy as np, gudhi
def torus_pts(n,seed,Rmaj=2.0,rmin=0.8):
    rng=np.random.default_rng(seed)
    th=rng.uniform(0,2*np.pi,n); ph=rng.uniform(0,2*np.pi,n)
    return np.column_stack([(Rmaj+rmin*np.cos(th))*np.cos(ph),(Rmaj+rmin*np.cos(th))*np.sin(ph),rmin*np.sin(th)])
def sphere_pts(n,seed):
    rng=np.random.default_rng(seed); v=rng.normal(size=(n,3)); v/=np.linalg.norm(v,axis=1,keepdims=True); return v
def top_bars(points,name):
    st=gudhi.AlphaComplex(points=points).create_simplex_tree(); st.compute_persistence()
    print(name)
    for d in (0,1,2):
        iv=st.persistence_intervals_in_dimension(d)
        lifes=sorted([( (1e9 if np.isinf(de) else de)-bi ) for (bi,de) in iv],reverse=True)
        ninf=sum(1 for (bi,de) in iv if np.isinf(de))
        # count "clearly persistent": lifetime > 10x the largest noise bar
        finite=sorted([de-bi for (bi,de) in iv if not np.isinf(de)],reverse=True)
        print(f"  dim {d}: #inf_bars={ninf}, top finite lifetimes={[round(x,3) for x in finite[:4]]}")
if __name__=="__main__":
    top_bars(sphere_pts(600,0),"S^2 (expect b0=1 inf, b1=0, b2=1 one long bar):")
    top_bars(torus_pts(900,0), "T^2 (expect b0=1 inf, b1=2 long bars, b2=1 long bar):")
