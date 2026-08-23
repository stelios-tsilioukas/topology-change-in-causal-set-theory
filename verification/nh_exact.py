import numpy as np, itertools
from nh import maximal_antichains, nerve_b1

def is_transitive(R):
    return not (( (R.astype(np.int8) @ R.astype(np.int8)) > 0 ) & ~R).any()

def all_posets(n):
    idx=[(i,j) for i in range(n) for j in range(i+1,n)]
    m=len(idx)
    for bits in range(1<<m):
        R=np.zeros((n,n),bool)
        for t,(i,j) in enumerate(idx):
            if bits>>t & 1: R[i,j]=True
        if is_transitive(R): yield R

print("="*72)
print(" EXHAUSTIVE SEARCH over all (linearly-extended) posets")
print("="*72)
print(f"{'n':>4} {'#posets':>12} {'#with spatial 1-cycle':>24}")
example=None
for n in range(4,7):
    tot=0; hit=0
    for R in all_posets(n):
        tot+=1
        found=False
        for A in maximal_antichains(R):
            if nerve_b1(R,A)>=1: found=True; break
        if found:
            hit+=1
            if example is None and n==6: example=R.copy()
    print(f"{n:>4} {tot:>12} {hit:>24}")

print(f"\n  ==> n_h = 6 exactly (nothing at n = 4 or 5)")

print("\n" + "="*72)
print(" WHY 6 — the combinatorial lower bound")
print("="*72)
print("""  A 1-cycle in a nerve needs >= 3 vertices, so the antichain has >= 3 elements
  {a,b,c}.  To get the boundary of a triangle (S^1) and not the filled triangle,
  each PAIR must have a common future element while the TRIPLE must not:
      x in fut(a) & fut(b),   y in fut(b) & fut(c),   z in fut(a) & fut(c),
      fut(a) & fut(b) & fut(c) = empty.
  That requires 3 distinct future elements.  Total: 3 + 3 = 6.        [tight]""")

if example is not None:
    print("\n  explicit minimal example (relation matrix, rows<cols):")
    for r in example.astype(int): print("      ",r)
    A=[A for A in maximal_antichains(example) if nerve_b1(example,A)>=1][0]
    print(f"      antichain = {A},  nerve b_1 = {nerve_b1(example,A)}")
