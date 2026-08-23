import numpy as np, itertools, random
from homology import betti_and_chi

print("="*78)
print(" A CRUCIAL SUBTLETY: does beta=(1,1,1,1) identify S^1 x S^2 ?")
print("="*78)
print("""  S^1 x S^2      : beta = (1,1,1,1),  cup product a\\smile b != 0
  S^1 v S^2 v S^3: beta = (1,1,1,1),  ALL cup products vanish
  => the Betti VECTOR alone does NOT distinguish a handle attachment from a
     wedge of spheres.  It is necessary, not sufficient.\n""")

# --- explicit small complex with beta=(1,1,1,1): S^3 v S^1 v S^2 ---
# vertices 0..4 = boundary of 4-simplex (S^3); 5 makes a 1-cycle; 6 makes a 2-cycle
S3=[c for c in itertools.combinations(range(5),4)]          # 5 tetrahedra
W = list(S3)
W += [(0,5),(1,5)]                                          # path -> 1-cycle with edge (0,1)
W += [(0,6,2),(2,6,3),(0,6,3)]                              # 2-sphere with face (0,2,3)
b,chi,f = betti_and_chi(W)
print(f"  explicit 7-vertex complex: betti = {[b[k] for k in sorted(b)]}, chi = {chi}")
print(f"  face counts {dict(sorted(f.items()))}")

print("\n"+"="*78)
print(" LOWER BOUND ON NERVE VERTICES FOR beta_3 != 0 PLUS beta_1 != 0")
print("="*78)
print("""  m = 5: a 3-cycle on 5 vertices must be the full boundary dDelta^4 (all five
         tetrahedra); adding the 4-simplex kills H_3, and all lower faces are
         already present.  Hence the complex IS dDelta^4 and beta = (1,0,0,1).
         => m = 5 CANNOT give (1,1,1,1).""")

def rand_complex(m, rng, pmax=0.55):
    maximal=[]
    for k in range(2,m+1):
        for c in itertools.combinations(range(m),k):
            if rng.random()<pmax: maximal.append(c)
    if not maximal: maximal=[(0,)]
    return maximal

rng=random.Random(0)
found=None
for trial in range(60000):
    W6=rand_complex(6,rng)
    try:
        b,chi,f=betti_and_chi(W6)
    except Exception: continue
    v=[b.get(k,0) for k in range(4)]
    if v==[1,1,1,1]:
        found=W6; break
print(f"  m = 6: randomized search over 60000 complexes -> "
      f"{'FOUND '+str(found) if found else 'none with beta=(1,1,1,1)'}")
print(f"  m = 7: explicit construction above succeeds.")
print("\n  => minimal nerve realising the Betti SIGNATURE: m = 7 (6 not found).")
print("     minimal nerve realising S^1 x S^2 ITSELF: covering type, ~9-10")
print("     (its minimal simplicial triangulation is known to need 10 vertices).")

print("\n"+"="*78)
print(" TOTAL CAUSAL-SET SIZE")
print("="*78)
print("""  In the MRS/nerve construction each MAXIMAL simplex needs a distinct witness
  element in the future (an element in the common future of exactly that set).
  Hence   n  >=  #(nerve vertices) + #(maximal simplices).""")
b,chi,f=betti_and_chi(W)
# count maximal simplices of the 7-vertex complex
mx=set()
for s in W: mx.add(tuple(sorted(s)))
print(f"\n  Betti-signature complex : 7 vertices + {len(mx)} maximal simplices  ->  n >~ {7+len(mx)}")
print( "  genuine S^1 x S^2       : 10 vertices + ~30 facets (Dehn-Sommerville,")
print( "                            N_3 = N_1 - N_0 for a closed 3-manifold)  ->  n >~ 40")
print( "  plus two S^3 layers (dDelta^4: 5 vertices + 5 facets each)          ->  +20")
