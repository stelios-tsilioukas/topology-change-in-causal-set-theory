import numpy as np
from itertools import combinations

def random_causet(n, p, seed):
    """Transitive percolation: relate i<j with prob p, then transitive closure."""
    rng=np.random.default_rng(seed)
    R=np.zeros((n,n),bool)
    for i in range(n):
        for j in range(i+1,n):
            if rng.random()<p: R[i,j]=True
    # transitive closure (Floyd-Warshall style)
    for k in range(n):
        R |= (R[:,[k]] & R[[k],:])
    return R

def chain_counts(R, kmax=6):
    """c[j] = number of j-element chains (totally ordered subsets)."""
    n=R.shape[0]
    c={1:n}
    # dynamic programming over chain length
    # paths[j][v] = number of j-chains ending at v
    paths={1:np.ones(n,dtype=object)}
    for j in range(2,kmax+1):
        prev=paths[j-1]
        cur=np.zeros(n,dtype=object)
        for v in range(n):
            preds=np.nonzero(R[:,v])[0]
            cur[v]=sum(prev[u] for u in preds)
        paths[j]=cur
        c[j]=int(cur.sum())
        if c[j]==0: break
    return c

def chi_order_complex(R,kmax=8):
    c=chain_counts(R,kmax)
    # order complex: (j)-element chain = (j-1)-simplex ; chi = sum_j (-1)^{j-1} c_j
    return sum((-1)**(j-1)*c[j] for j in sorted(c))

print("=== Is the order-complex Euler characteristic order-reversal invariant? ===\n")
print(f"{'n':>4} {'p':>5} {'chi(C)':>12} {'chi(C*)':>12}  equal?")
for (n,p,seed) in [(12,0.3,0),(12,0.5,1),(16,0.25,2),(16,0.4,3),(20,0.2,4),(20,0.35,5),(24,0.3,6)]:
    R=random_causet(n,p,seed)
    Rd=R.T.copy()          # order reversal: x <* y  iff  y < x
    a=chi_order_complex(R); b=chi_order_complex(Rd)
    print(f"{n:>4} {p:>5} {a:>12} {b:>12}   {a==b}")

print("\n=== Why: the chain sets are identical ===")
R=random_causet(14,0.35,7); Rd=R.T.copy()
ca=chain_counts(R); cb=chain_counts(Rd)
print(" j-chain counts  C :", {k:v for k,v in ca.items()})
print(" j-chain counts  C*:", {k:v for k,v in cb.items()})
print(" identical:", ca==cb)

print("\n=== Morse-index bookkeeping in n=4 ===")
for lam in (0,1,2,3,4):
    print(f"  index {lam}: (-1)^{lam} = {(-1)**lam:+d}   ->  time-reversed index {4-lam}: (-1)^{4-lam} = {(-1)**(4-lam):+d}")
print("\n  => in EVEN dimension n=4, (-1)^lambda is INVARIANT under lambda -> n-lambda.")
print("  (in odd n it would flip sign)")

print("\n=== Handle decomposition of the two species ===")
print("  Euclidean wormhole  M # (S1xS3): dchi = 0 - 2 = -2  =  (-1)^1 + (-1)^3 = -1 -1")
print("     -> ONE index-1 and ONE index-3 critical point")
print("  Nariai instanton    M # (S2xS2): dchi = 4 - 2 = +2  =  (-1)^2 + (-1)^2 = +1 +1")
print("     -> TWO index-2 critical points")
print("\n  Time reversal maps index 1 <-> 3 (wormhole -> wormhole) and 2 <-> 2 (Nariai -> Nariai).")
print("  => each species maps to ITSELF. No pairing of +2 with -2 exists.")
