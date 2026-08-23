import numpy as np
from fractions import Fraction

def chain_counts(R, kmax=40):
    """c[j] = number of j-element chains in the poset R (strict order)."""
    n=R.shape[0]
    paths={1:[1]*n}          # paths[j][v] = # j-chains ending at v
    c={1:n}
    preds=[np.nonzero(R[:,v])[0] for v in range(n)]
    for j in range(2,kmax+1):
        prev=paths[j-1]
        cur=[sum(prev[u] for u in preds[v]) for v in range(n)]
        tot=sum(cur)
        if tot==0: break
        paths[j]=cur; c[j]=tot
    return c

def chi_order_complex(R):
    """chi(Delta(P)) = sum_j (-1)^{j-1} c_j   ((j)-chain = (j-1)-simplex)"""
    c=chain_counts(R)
    return sum((-1)**(j-1)*c[j] for j in sorted(c))

# ---------- geometries ----------
def diamond2d(N,seed,with_tips=True):
    rng=np.random.default_rng(seed)
    # causal diamond in 1+1D: light-cone coords u,v in (0,1)
    u=rng.random(N); v=rng.random(N)
    if with_tips:
        u=np.concatenate(([0.0,1.0],u)); v=np.concatenate(([0.0,1.0],v))
    n=len(u)
    R=(u[:,None]<u[None,:])&(v[:,None]<v[None,:])
    np.fill_diagonal(R,False)
    return R

def cylinder2d(N,seed,L=1.0,T=0.30):
    """1+1D cylinder slab: spatial circle circumference L, time extent T.
       manifold = S^1 x I  ->  chi = 0."""
    rng=np.random.default_rng(seed)
    t=rng.random(N)*T; x=rng.random(N)*L
    dt=t[None,:]-t[:,None]
    dx=np.abs(x[None,:]-x[:,None]); dx=np.minimum(dx,L-dx)
    R=(dt>0)&(dt>dx)
    np.fill_diagonal(R,False)
    return R

print("="*72)
print(" DOES THE ORDER-COMPLEX EULER CHARACTERISTIC RECOVER THE MANIFOLD chi?")
print("="*72)

print("\n(A) 1+1D causal diamond WITH tips  (manifold: 2-ball, chi = +1)")
for N in (20,40,80,160):
    vals=[chi_order_complex(diamond2d(N,s,True)) for s in range(4)]
    print(f"    N={N:>4}:  chi_oc = {vals}")

print("\n(B) 1+1D causal diamond WITHOUT tips  (sprinkling only)")
for N in (20,40,80,160):
    vals=[chi_order_complex(diamond2d(N,s,False)) for s in range(4)]
    print(f"    N={N:>4}:  chi_oc = {[f'{v:+d}' for v in vals]}")

print("\n(C) 1+1D cylinder slab S^1 x I  (manifold chi = 0)")
for N in (40,80,160,320):
    vals=[chi_order_complex(cylinder2d(N,s)) for s in range(4)]
    print(f"    N={N:>4}:  chi_oc = {[f'{v:+d}' for v in vals]}")
