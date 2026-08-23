import numpy as np
print("="*76)
print(" 1. HANDLE DECOMPOSITION CHECK  (is delta_chi = -2 consistent?)")
print("="*76)
print("""  S^1 x S^3 has handle decomposition: 0-handle, 1-handle, 3-handle, 4-handle
      chi = (+1) + (-1) + (-1) + (+1) = 0   [= chi(S^1)chi(S^3) = 0]  OK

  Connected sum M # (S^1 x S^3) absorbs the 0- and 4-handles, leaving M with
  a 1-handle and a 3-handle attached:
      delta_chi = (-1)^1 + (-1)^3 = -2      [= chi(M1)+chi(M2)-2]      OK

  Cobordism W: S^3 -> S^3 containing that pair:
      chi(W) = chi(S^3) + (-1)^1 + (-1)^3 = 0 - 2 = -2
  and W = (S^1 x S^3) minus two 4-balls:
      chi = 0 - 2*(1) = -2                                             OK
  Capping the two S^3 boundaries with 4-balls returns chi = -2 + 2 = 0. OK
  => the -2 (handle pair) and the 0 (closed manifold) are consistent.""")

print("\n"+"="*76)
print(" 2. SPATIAL BETTI VECTORS")
print("="*76)
S3=(1,0,0,1); S1S2=(1,1,1,1)
chi=lambda b: sum((-1)**k*x for k,x in enumerate(b))
print(f"   S^3        beta = {S3}   chi = {chi(S3)}")
print(f"   S^1 x S^2  beta = {S1S2}   chi = {chi(S1S2)}")
print(f"   Delta beta (creation)    = {tuple(a-b for a,b in zip(S1S2,S3))}")
print(f"   Delta chi  (spatial)     = {chi(S1S2)-chi(S3)}   <-- carries NO information")
print("   => the Betti VECTOR is the discriminating spatial observable, not chi.")

print("\n"+"="*76)
print(" 3. CONSEQUENCE FOR THE SELECTED SCALE  (the step they did not take)")
print("="*76)
p=1e-4; pref=32*np.pi**2; Lam=1e-122
def L_of(nh): return (pref*p/Lam)**(1.0/(4*nh))
print(f"{'n_h':>6} {'L = l_cg/l_P':>15} {'E = M_Pl/L [GeV]':>20} {'Higgs bound (L>=37)':>22}")
for nh in (6,10,15,19,20,25,30,45):
    L=L_of(nh); E=1.22e19/L
    ok = "OK" if L>=37 else "VIOLATED"
    print(f"{nh:>6} {L:>15.3g} {E:>20.2e} {ok:>22}")
from scipy.optimize import brentq
nstar=brentq(lambda n: L_of(n)-37.0, 2, 200)
print(f"\n   Higgs triviality bound (l_k >= 37 l_P) requires   n_h <= {nstar:.1f}")
print(f"   Covering-type argument for S^1 x S^2 suggests      n_h >~ 19-20")
print(f"   Explicit star-cover realisation gives              n_h <~ 45")
print("\n   => the window is NARROW and the upper part of the range is EXCLUDED.")
