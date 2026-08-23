import numpy as np
from scipy.optimize import brentq
Mpl=1.22e19; Mp=Mpl/np.sqrt(8*np.pi); coef=np.sqrt(6)*np.pi/2
A0=8*np.pi/(2*np.sqrt(6)*np.pi**2*coef)      # a_0^2/l_P^2 = A0 * S
pref=8.0/A0**2                                # Lambda l_P^2 = pref * alphabar * e^-S
print("="*78); print(" alpha-DEPENDENCE OF THE PREDICTION  (alphabar = alpha/l_P^2)"); print("="*78)
print(f"   a_0^2/l_P^2 = {A0:.4f} S        Lambda_eff l_P^2 = {pref:.1f} * alphabar * e^-S")
print(f"   =>  S_1 = ln({pref:.1f} alphabar) + 281.0 = {np.log(pref)+280.9:.1f} + ln(alphabar)\n")
print(f"{'alphabar':>12} {'S_1':>9} {'f [GeV]':>13} {'a_0/l_P':>10}")
for ab in (1e-4,1e-2,1.0,1e2,1e4):
    S=brentq(lambda S: np.log(pref*ab)-S-np.log(1e-122),1,2000)
    f=coef*Mp/S; a0=np.sqrt(A0*S)
    print(f"{ab:>12.0e} {S:>9.1f} {f:>13.3e} {a0:>10.2f}")
print("""
   Varying the Gauss-Bonnet coupling over EIGHT orders of magnitude changes the
   predicted decay constant by ~6 per cent: alpha enters only logarithmically,
   because it multiplies the prefactor while S sits in an exponent.""")
print("\n"+"="*78); print(" ROUND-HANDLE (resolution-selected) CASE, for comparison"); print("="*78)
c=3*np.pi**2/2
print(f"{'alphabar':>12} {'L = l_cg/l_P':>14} {'E [GeV]':>13}")
for ab in (1e-2,1.0,1e2):
    L=brentq(lambda L: np.log(32*np.pi**2*ab)-4*np.log(L)-c*L**2-np.log(1e-122),1.0001,1e6)
    print(f"{ab:>12.0e} {L:>14.2f} {1.22e19/L:>13.2e}")
