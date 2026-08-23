import numpy as np
from itertools import combinations
from fcnerve import diamond, cylinder

def masks(R):
    n=R.shape[0]
    Fut=[1<<x for x in range(n)]; Pas=[1<<x for x in range(n)]
    for x in range(n):
        for y in np.nonzero(R[x,:])[0]: Fut[x]|=1<<int(y)
        for y in np.nonzero(R[:,x])[0]: Pas[x]|=1<<int(y)
    return Fut,Pas

def interval_complex_chi(R, ncg):
    """Cover by ALEXANDROV INTERVALS of cardinality <= ncg (bounded regions!).
       S is a face iff S is contained in some such interval."""
    n=R.shape[0]; Fut,Pas=masks(R)
    faces=set()
    # all small intervals (including degenerate ones = single elements & links)
    covers=[]
    for x in range(n):
        covers.append(1<<x)
        for y in np.nonzero(R[x,:])[0]:
            M=Fut[x]&Pas[int(y)]
            if bin(M).count("1")<=ncg: covers.append(M)
    covers=set(covers)
    # mark every subset of every cover element
    for M in covers:
        bits=[i for i in range(n) if M>>i & 1]
        for k in range(1,len(bits)+1):
            for c in combinations(bits,k):
                s=0
                for b in c: s|=1<<b
                faces.add(s)
    chi=0; dims={}
    for S in faces:
        k=bin(S).count("1"); chi+=(-1)**(k-1); dims[k]=dims.get(k,0)+1
    return chi,dims,len(covers)

print("="*74)
print(" NERVE OF A COVER BY BOUNDED ALEXANDROV INTERVALS (|I| <= n_cg)")
print("="*74)
print(" The scale n_cg is FORCED: unbounded future cones give a trivial nerve.\n")

for ncg in (4,6,8):
    print(f" --- n_cg = {ncg} ---")
    print(f"   (A) diamond+tips  [chi=+1]: ", end="")
    print([interval_complex_chi(diamond(N,s),ncg)[0] for N,s in [(14,0),(14,1),(18,2),(18,3)]])
    print(f"   (B) cylinder slab [chi= 0]: ", end="")
    print([interval_complex_chi(cylinder(N,s),ncg)[0] for N,s in [(14,0),(14,1),(18,2),(18,3),(20,4)]])
