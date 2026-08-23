import numpy as np, itertools, sys
from homology import betti_and_chi

def close(R):
    n=R.shape[0]
    for k in range(n): R|=np.outer(R[:,k],R[k,:])
    return np.triu(R,1)

def rand_causet(n,p,rng):
    return close(np.triu(rng.random((n,n))<p,1))

def maximal_antichains(R):
    """all maximal antichains (Bron-Kerbosch on the incomparability graph)."""
    n=R.shape[0]
    comp = R | R.T
    adj=[set(np.nonzero(~comp[i])[0])-{i} for i in range(n)]
    out=[]
    def bk(Rset,P,X):
        if not P and not X:
            if len(Rset)>=3: out.append(sorted(Rset))
            return
        if not P: return
        piv=max(P|X, key=lambda v: len(adj[v]&P))
        for v in list(P-adj[piv]):
            bk(Rset|{v}, P&adj[v], X&adj[v])
            P=P-{v}; X=X|{v}
    bk(set(),set(range(n)),set())
    return out

def nerve_b1(R, A):
    """Cech nerve of future-sets of antichain elements (within the causet).
       simplex {a_i...} iff their strict futures intersect."""
    fut={a:set(np.nonzero(R[a,:])[0]) for a in A}
    m=len(A)
    if m<3: return 0
    maximal=[]
    # build all simplices whose futures have common intersection
    for k in range(m,0,-1):
        for c in itertools.combinations(range(m),k):
            inter=set.intersection(*[fut[A[i]] for i in c]) if k>1 else fut[A[c[0]]]
            if inter:
                # keep only if not already contained in a recorded maximal simplex
                if not any(set(c)<=set(s) for s in maximal):
                    maximal.append(c)
    cov=set()
    for s in maximal: cov.update(s)
    for i in range(m):
        if i not in cov: maximal.append((i,))
    if not maximal: return 0
    b,chi,f=betti_and_chi(maximal)
    return b.get(1,0)

def scan(n, trials, rng, ps=(0.25,0.35,0.45)):
    hits=0
    for t in range(trials):
        p=ps[t%len(ps)]
        R=rand_causet(n,p,rng)
        for A in maximal_antichains(R):
            if nerve_b1(R,A)>=1:
                hits+=1; break
    return hits/trials

rng=np.random.default_rng(7)
print("="*72)
print(" MINIMAL TOPOLOGY-CHANGING SUB-ORDER:  spatial 1-cycle in the MRS nerve")
print("="*72)
print("\n criterion: exists a maximal antichain whose future-set nerve has b_1 >= 1")
print(" (a spatial loop -> the handle is born; being finite, it must also die)\n")
print(f"{'n':>4} {'fraction with b1>=1':>22}")
TR=4000
first=None
for n in range(4,13):
    fr=scan(n,TR,rng)
    print(f"{n:>4} {fr:>22.4f}")
    if fr>0 and first is None: first=n
    sys.stdout.flush()
print(f"\n  minimal n with a spatial 1-cycle:  n_h = {first}")
