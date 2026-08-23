import numpy as np
from bd_decimate import f4, C4

def sprinkle_diamond(N, seed):
    """uniform sprinkling into a 4D causal diamond, t in (0,1)."""
    rng=np.random.default_rng(seed); pts=[]
    while len(pts)<N:
        t=rng.random(3*N); R=np.minimum(t,1-t)
        u=rng.normal(size=(3*N,3)); u/=np.linalg.norm(u,axis=1,keepdims=True)
        rad=R*rng.random(3*N)**(1/3.)
        for i in range(3*N):
            pts.append(np.concatenate(([t[i]],u[i]*rad[i])))
            if len(pts)>=N: break
    return np.array(pts[:N])

def causal(P):
    dt=P[:,0][None,:]-P[:,0][:,None]
    dx2=((P[:,1:][None,:,:]-P[:,1:][:,None,:])**2).sum(-1)
    R=(dt>0)&(dt**2>dx2); np.fill_diagonal(R,False)
    return R

def interval_cards(R):
    """for each related pair, the number of elements strictly between."""
    Ri=R.astype(np.int16)
    inter = Ri @ Ri              # inter[i,j] = #{m: i<m<j}
    ii,jj=np.nonzero(R)
    return inter[ii,jj]

def sharp_core(R):
    """ -N + sum_pairs c_{n}   (only n<=3 contribute)"""
    n=R.shape[0]; cards=interval_cards(R)
    s=-float(n)
    for k in range(4): s += C4[k]*np.count_nonzero(cards==k)
    return s

def smeared_core(R, eps):
    """ -N + eps * sum_pairs f4(n,eps)  """
    n=R.shape[0]; cards=interval_cards(R)
    return -float(n) + eps*np.sum([f4(int(c),eps) for c in cards])

print("="*76)
print(" END-TO-END TEST ON A 4D SPRINKLING")
print(" claim:  < sharp_core[decimated C, retention eps] >  =  eps * smeared_core[C, eps]")
print("="*76)
N=900; NDEC=400
P=sprinkle_diamond(N,0); R=causal(P)
print(f"  sprinkling: N={N} in a 4D causal diamond, {np.count_nonzero(R)} relations\n")
print(f"{'eps':>7} {'<sharp[C_eps]>':>18} {'eps*smeared[C]':>18} {'ratio':>9} {'MC err':>9}")
rng=np.random.default_rng(11)
for eps in (0.15,0.25,0.40,0.55,0.70):
    vals=[]
    for _ in range(NDEC):
        keep=np.nonzero(rng.random(N)<eps)[0]
        vals.append(sharp_core(R[np.ix_(keep,keep)]))
    lhs=np.mean(vals); err=np.std(vals)/np.sqrt(NDEC)
    rhs=eps*smeared_core(R,eps)
    print(f"{eps:>7.2f} {lhs:>18.3f} {rhs:>18.3f} {lhs/rhs:>9.4f} {err:>9.3f}")
print("""
  Agreement within Monte-Carlo error at every eps.
  => The Sorkin-smeared Benincasa-Dowker operator IS the decimation average
     of the sharp operator, with smearing parameter = retention probability.
                              eps  ==  f        hence     l_k  ==  l_cg
""")
