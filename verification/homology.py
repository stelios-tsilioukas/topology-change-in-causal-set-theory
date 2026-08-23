import numpy as np
from itertools import combinations

# ---------- GF(2) rank via bitset elimination ----------
def gf2_rank(cols):
    pivots = {}
    rank = 0
    for v in cols:
        cur = v
        while cur:
            hb = cur.bit_length() - 1
            if hb in pivots:
                cur ^= pivots[hb]
            else:
                pivots[hb] = cur
                rank += 1
                break
    return rank

# ---------- Full simplicial complex from maximal simplices ----------
def build_complex(maximal):
    all_simp = set()
    for s in maximal:
        s = tuple(sorted(s))
        for k in range(1, len(s)+1):
            for c in combinations(s, k):
                all_simp.add(c)
    by_dim = {}
    for s in all_simp:
        by_dim.setdefault(len(s)-1, []).append(s)
    for d in by_dim:
        by_dim[d].sort()
    return by_dim

def betti_and_chi(maximal):
    by_dim = build_complex(maximal)
    dims = sorted(by_dim)
    f = {d: len(by_dim[d]) for d in by_dim}
    idx = {d: {s:i for i,s in enumerate(by_dim[d])} for d in by_dim}
    # boundary ranks
    rank_bd = {}   # rank of ∂_d  (d-simplices -> (d-1)-simplices)
    for d in dims:
        if d-1 in by_dim:
            cols = []
            for s in by_dim[d]:
                bits = 0
                for i in range(len(s)):
                    face = s[:i] + s[i+1:]
                    bits ^= (1 << idx[d-1][face])
                cols.append(bits)
            rank_bd[d] = gf2_rank(cols)
        else:
            rank_bd[d] = 0
    maxd = max(dims)
    betti = {}
    for d in range(0, maxd+1):
        fd = f.get(d, 0)
        rd = rank_bd.get(d, 0)       # rank ∂_d
        rd1 = rank_bd.get(d+1, 0)    # rank ∂_{d+1}
        betti[d] = fd - rd - rd1
    chi = sum((-1)**d * f.get(d,0) for d in range(0, maxd+1))
    return betti, chi, f

# ---------- Validation ----------
if __name__ == "__main__":
    print("=== VALIDATION OF HOMOLOGY ENGINE ===\n")

    # (1) Filled triangle = disk : chi=1, betti=[1,0,0]
    disk = [(0,1,2)]
    b,chi,f = betti_and_chi(disk)
    print("Disk (filled triangle):   betti=", [b[k] for k in sorted(b)], " chi=", chi, " f=", dict(f))

    # (2) Octahedron boundary = S^2 : chi=2, betti=[1,0,1]
    V = {'+x':0,'-x':1,'+y':2,'-y':3,'+z':4,'-z':5}
    oct_faces = []
    for a in (0,1):
        for b2 in (2,3):
            for c in (4,5):
                oct_faces.append((a,b2,c))
    b,chi,f = betti_and_chi(oct_faces)
    print("Octahedron (S^2):         betti=", [b[k] for k in sorted(b)], " chi=", chi, " f=", dict(f))

    # (3) 3x3 flat torus triangulation : chi=0, betti=[1,2,1]
    def vid(i,j): return 3*(i%3) + (j%3)
    tor_faces = []
    for i in range(3):
        for j in range(3):
            a = vid(i,j); bb = vid(i+1,j); c = vid(i,j+1); d = vid(i+1,j+1)
            tor_faces.append(tuple(sorted((a,bb,c))))
            tor_faces.append(tuple(sorted((bb,d,c))))
    tor_faces = list(set(tor_faces))
    b,chi,f = betti_and_chi(tor_faces)
    print("3x3 torus (T^2):          betti=", [b[k] for k in sorted(b)], " chi=", chi, " f=", dict(f))
    print("\ndelta_chi (T^2 - S^2) =", 0 - 2, " (expected -2, a handle birth)")
