import numpy as np
from scipy.optimize import brentq

print("="*78)
print(" FOAM SIZE SPECTRUM FROM THE PATH INTEGRAL")
print("="*78)
print("""
 A handle of physical size l_w has Euclidean action, by dimensions,

     S(l_w) = c (l_w/l_P)^2          [ ~ (1/16 pi G) * R * Vol ~ a^2/G ]

 the Gauss-Bonnet piece being topological and hence size-independent.
 Nucleation weight  ~ e^{-S}, positions ~ V/l_w^4, so

     n(l_w) d ln l_w  ~  l_w^{-4} e^{-c (l_w/l_P)^2} d ln l_w .

 At coarse-graining scale l_cg only handles with l_w >~ l_cg are RESOLVED
 (established numerically: detection is a size threshold, not a pattern count).
 The steeply falling spectrum makes the integral dominated by l_w ~ l_cg:

     dchi/dV|_cg  ~  (dchi / l_cg^4) e^{-S(l_cg)}

 and hence, with alpha = l_P^2 and |dchi| = 2,
""")
print("     Lambda_eff l_P^2  =  32 pi^2 L^{-4} e^{-c L^2},     L = l_cg/l_P\n")

Lam=1e-122; pref=32*np.pi**2
def eq(L,c): return np.log(pref)-4*np.log(L)-c*L**2-np.log(Lam)
print("="*78); print(" SOLVING FOR THE SELECTED SCALE"); print("="*78)
print(f"{'c':>8} {'L = l_cg/l_P':>14} {'E = M_Pl/L [GeV]':>20} {'S(l_cg)':>10} {'Higgs L>=37':>13}")
for c in (0.01,0.02,1/(16*np.pi),0.1,0.199,0.4,1.0):
    L=brentq(lambda x: eq(x,c), 1.0001, 1e6)
    print(f"{c:>8.4f} {L:>14.1f} {1.22e19/L:>20.2e} {c*L**2:>10.1f} "
          f"{'OK' if L>=37 else 'violated':>13}")
ccrit=brentq(lambda c: brentq(lambda x: eq(x,c),1.0001,1e6)-37.0, 1e-4, 5)
print(f"\n  Higgs bound (L >= 37) requires   c <= {ccrit:.3f}")

print("\n"+"="*78); print(" THE STRUCTURAL POINT"); print("="*78)
S_req = np.log(pref)-np.log(Lam)
print(f"""
 Setting Lambda_eff = Lambda_obs requires

     S(l_cg) = ln(32 pi^2 / Lambda_obs l_P^2) - 4 ln L  ~  {S_req:.0f} - 4 ln L  ~  280 .

 The 122 orders of magnitude are therefore NOT produced by a large power of a
 ratio, but by an exponent of order  ln(10^122) = {np.log(10**122 if False else 10)*122:.0f}.
 The statement becomes: the smallest RESOLVABLE handle has Euclidean action
 ~ 280 in Planck units.  This is the same structure as instanton suppression
 or dimensional transmutation, and is far less fine-tuned than a 10^-122 ratio.
""")
print("="*78); print(" CONSISTENCY OF THE REVIVED MECHANISM"); print("="*78)
print("""  * w = -1 exactly: p_h and l_cg are properties of the local order and of the
    action, neither epoch dependent.  (Sec. VI.C unchanged.)
  * l_cg = l_k still holds: the binomial thinning theorem is independent of
    the flow and survives intact.
  * The fluctuation channel remains excluded (independent argument).
  * Coleman's rate is RECOVERED, not replaced: the CST content is that the
    coarse-graining scale selects WHICH handle size dominates.""")
