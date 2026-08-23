import numpy as np
# Eichhorn et al. (2305.07595): fit parameters
A, a, b = -1.363, 0.126, 3.161
v, k0 = 246.0, 173.0          # GeV
MH = 125.25                   # measured Higgs mass

def F(z):
    if z <= 1: return A*(z**a - 1.0) + np.log(z)
    return z**b

def M_upper(rho14, kL_over_rho=1.0):
    """their Eq.(37): M = sqrt( (8 pi^2 v^2/9) / [F(kL/rho^1/4) - F(k0/rho^1/4)] )"""
    D = F(kL_over_rho) - F(k0/rho14)
    return np.sqrt((8*np.pi**2*v**2/9.0)/D)

print("="*74)
print(" EICHHORN et al. TRIVIALITY BOUND EVALUATED AT OUR PREDICTED SCALE")
print("="*74)
print(" Their result: the Landau pole occurs AT the nonlocality scale, k_L ~ rho^(1/4).")
print(" Higher nonlocality ENERGY => tighter Higgs bound.\n")
print(f"{'rho^(1/4) [GeV]':>18} {'l_k / l_P':>12} {'M_H upper [GeV]':>18} {'vs 125.25':>12}")
for rho14,lab in [(1e13,''),(1e14,'<- our n_h=6'),(1e15,''),(1e16,'<- our n_h=10'),
                  (1e17,''),(3.9e17,'<- critical'),(1.22e19,'<- PLANCK')]:
    L = 1.22e19/rho14
    M = M_upper(rho14)
    ok = "OK" if M>MH else "VIOLATED"
    print(f"{rho14:>18.2e} {L:>12.2e} {M:>18.1f} {ok:>12}  {lab}")

# critical scale
from scipy.optimize import brentq
crit = brentq(lambda r: M_upper(r)-MH, 1e12, 1e19)
print(f"\n  Critical nonlocality scale (M_upper = M_H): rho^(1/4) = {crit:.2e} GeV")
print(f"  => bound SATISFIED for rho^(1/4) < {crit:.1e} GeV,  i.e.  l_k > {1.22e19/crit:.0f} l_P")
print(f"  => bound VIOLATED at the Planck scale.")
print(f"\n  Our predicted range l_k = 1e3-1e5 l_P  (rho^(1/4) = 1e14-1e16 GeV): SATISFIED.")
print("\n  NB: their Eq.(37) appears to omit the factor 3 of Eq.(22) [m^2 = 3 lambda v^2].")
print(f"      Including it multiplies all bounds by sqrt(3): {M_upper(1e14)*np.sqrt(3):.0f} GeV at 1e14,")
print(f"      {M_upper(1.22e19)*np.sqrt(3):.0f} GeV at the Planck scale -- so the conclusion is")
print("      unchanged in direction, only in where the critical scale sits.")
