import numpy as np

def future_masks(R):
    """inclusive future bitmask of each element"""
    n=R.shape[0]
    F=[]
    for x in range(n):
        m=1<<x
        for y in np.nonzero(R[x,:])[0]: m|=1<<int(y)
        F.append(m)
    return F

def nerve_chi_and_counts(R, nmax=22):
    """chi of the nerve of inclusive future cones:
       S is a simplex iff  AND_{x in S} J+(x)  != 0 .
       chi = sum_{S != empty, face} (-1)^{|S|-1}"""
    n=R.shape[0]
    assert n<=nmax, "too big for exhaustive subset DP"
    F=future_masks(R)
    G=[0]*(1<<n); G[0]=(1<<n)-1
    chi=0; counts={}
    for S in range(1,1<<n):
        low=(S & -S).bit_length()-1
        G[S]=G[S ^ (1<<low)] & F[low]
        if G[S]:
            k=bin(S).count("1")
            chi += (-1)**(k-1)
            counts[k]=counts.get(k,0)+1
    return chi, counts

# ---------------- geometries (1+1D) ----------------
def diamond(N,seed,tips=True):
    rng=np.random.default_rng(seed)
    u=rng.random(N); v=rng.random(N)
    if tips: u=np.concatenate(([0.,1.],u)); v=np.concatenate(([0.,1.],v))
    R=(u[:,None]<u[None,:])&(v[:,None]<v[None,:]); np.fill_diagonal(R,False)
    return R

def cylinder(N,seed,L=1.0,T=0.45):
    rng=np.random.default_rng(seed)
    t=rng.random(N)*T; x=rng.random(N)*L
    dt=t[None,:]-t[:,None]
    dx=np.abs(x[None,:]-x[:,None]); dx=np.minimum(dx,L-dx)
    R=(dt>0)&(dt>dx); np.fill_diagonal(R,False)
    return R

print("="*74)
print(" NERVE OF INCLUSIVE FUTURE CONES  -- a purely order-theoretic 4D complex")
print("="*74)
print(" simplex {x_1..x_k}  <=>  the x_i have a common future element\n")

print(" (A) 1+1D causal diamond with tips   [manifold: 2-ball, chi = +1]")
for N in (8,12,16,18):
    v=[nerve_chi_and_counts(diamond(N,s))[0] for s in range(4)]
    print(f"     N={N+2:>3}:  chi_nerve = {v}")

print("\n (B) 1+1D cylinder slab S^1 x I     [manifold chi = 0, homotopy S^1]")
for N in (10,14,18,20):
    v=[nerve_chi_and_counts(cylinder(N,s))[0] for s in range(5)]
    print(f"     N={N:>3}:  chi_nerve = {v}")
