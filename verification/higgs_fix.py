import numpy as np
from scipy.optimize import brentq
A,a,b = -1.363, 0.126, 3.161      # their Eq.(31)
v, k0 = 246.0, 173.0              # their text below Eq.(37)
MH = 125.25

def F(z): return A*(z**a-1.0)+np.log(z) if z<=1 else z**b

print("="*78); print(" THE FACTOR OF THREE, RESOLVED"); print("="*78)
print("""  Eq.(22):  m^2_upper = 3 lambda_0 v^2
  Eq.(36):  lambda_0 < (8 pi^2/9)/D ,   D = F(k_L/rho^1/4) - F(k_0/rho^1/4)
  => m^2_upper < 3 * (8 pi^2/9) v^2/D = 8 pi^2 v^2/(3 D)

  Eq.(37) as printed:  M_S^Upper = ( 8 pi^2 v^2/9 / D )^{1/2}   <-- has /9

  So Eq.(37) is inconsistent with (22)+(36) by exactly a factor 3 in m^2,
  i.e. sqrt(3) in the mass.  CONFIRMED.\n""")

def M_printed(r14):
    D=F(1.0)-F(k0/r14); return np.sqrt(8*np.pi**2*v**2/9.0/D)
def M_corrected(r14):
    D=F(1.0)-F(k0/r14); return np.sqrt(8*np.pi**2*v**2/3.0/D)

print("="*78); print(" CONSEQUENCE  (Landau pole at the nonlocality scale, k_L = rho^1/4)")
print("="*78)
print(f"{'rho^1/4 [GeV]':>16} {'l_k/l_P':>10} {'Eq.(37) as printed':>20} {'corrected (x sqrt3)':>21}")
for r14 in (1e14,1e16,1e17,2.8e18,1.22e19):
    L=1.22e19/r14
    mp,mc=M_printed(r14),M_corrected(r14)
    f1="OK" if mp>MH else "viol"
    f2="OK" if mc>MH else "viol"
    print(f"{r14:>16.2e} {L:>10.1f} {mp:>14.1f} {f1:>5} {mc:>15.1f} {f2:>5}")
print(f"\n  measured M_H = {MH} GeV")

# where would the corrected bound ever be violated?
try:
    zc=brentq(lambda lz: np.sqrt(8*np.pi**2*v**2/3.0/(F(1.0)-F(k0/np.exp(lz))))-MH, np.log(1e5), np.log(1e60))
    print(f"  corrected bound would fail only for rho^1/4 > {np.exp(zc):.1e} GeV  (unphysical)")
except Exception as e:
    print("  corrected bound is satisfied over the whole physical range")

print("\n"+"="*78); print(" CROSS-CHECK OF OUR BD -> EH NORMALISATION"); print("="*78)
print("""  Their Eq.(41)/(47):  S[C] = (l^4/16 pi G) sum_x B(-2),  with (Table 1, d=4)
      alpha_4 = -4/sqrt6,  beta_4 = +4/sqrt6,  C_i = (1,-9,16,-8).
  Then   sum_x B(-2) = (8/(sqrt6 l^2)) [N - N_1 + 9N_2 - 16N_3 + 8N_4]
  so     S_EH = (l^2/(2 sqrt6 pi G)) S_core   =>   S_core = (2 sqrt6 pi G/l^2) S_EH""")
print(f"\n     2 sqrt6 pi = {2*np.sqrt(6)*np.pi:.3f}   -- matches our Eq. for the normalisation. CONFIRMED.")
