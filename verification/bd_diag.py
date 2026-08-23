import numpy as np
from bd_handle import sprinkle_box, causal_flat, bd_core_smeared

print("="*76)
print(" DIAGNOSIS: why the direct BD handle calculation is out of reach")
print("="*76)
print("\n (a) flat-space baseline (R=0, so the BULK action must vanish):")
for N in (300,500,700,900):
    v=[]
    for s in range(6):
        t,x=sprinkle_box(N,2.0,2.0,s)
        v.append(bd_core_smeared(causal_flat(t,x),0.35))
    print(f"     N={N:>4}:  S_core = {np.mean(v):>9.1f} +/- {np.std(v):>6.1f}")
print("""     The value is large and grows with N: these are BOUNDARY terms of the
     finite sprinkling region, not bulk curvature.  A handle signal must be
     extracted on top of this, and it is far smaller.""")

print("\n (b) resolution requirement:")
print("""     A continuum handle needs  a >> l  (to be a geometry, not a few points)
     AND  box >> a  (to be a local perturbation).  Taking a = 3l and box = 5a
     gives box = 15 l, hence
          N = (box/l)^4 = 15^4 = 50625 elements,
     with O(N^2) relation storage (2.5x10^9) and O(N^2) interval counting.
     That is beyond this environment; the runs above have a/l <= 1, i.e. the
     'handle' is sub-discreteness and only its global causal reconnection is
     seen -- which is why the measured dS was ~constant in a.""")

print("\n"+"="*76)
print(" WHAT CAN BE DONE EXACTLY: the BD -> Einstein-Hilbert normalisation")
print("="*76)
print("""   From <B.1>(x) -> -R(x)/2  and  sum_x -> (1/l^4) int sqrt(g) d^4x :

        (4/(sqrt6 l^2)) S_core  =  -(1/2 l^4) int sqrt(g) R
        S_core = -(sqrt6/(8 l^2)) int sqrt(g) R

   and with the Euclidean  S_EH = -(1/16 pi G) int sqrt(g) R :""")
conv=2*np.sqrt(6)*np.pi
print(f"        S_core = (2 sqrt6 pi G / l^2) S_EH   ->   S_core = {conv:.3f} S_EH   (l^2 = G)")
print(f"        i.e.   S_EH = S_core / {conv:.3f}")

print("\n"+"="*76)
print(" CONTINUUM VALUE OF c FOR THE ROUND S^1 x S^3 WORMHOLE")
print("="*76)
print("""   ds^2 = dtau^2 + a^2 dOmega_3^2 ,  S^1 circumference 2 pi a
        R = 6/a^2          (S^3 factor;  S^1 is flat)
        Vol = (2 pi a)(2 pi^2 a^3) = 4 pi^3 a^4
        int sqrt(g) R = 24 pi^3 a^2""")
c=24*np.pi**3/(16*np.pi)
print(f"        S_EH = -(1/16 pi G) * 24 pi^3 a^2 = -(3 pi^2/2)(a^2/G)")
print(f"        |c| = 3 pi^2 / 2 = {3*np.pi**2/2:.2f}")

print("\n"+"="*76)
print(" CONSEQUENCE")
print("="*76)
from scipy.optimize import brentq
Lam=1e-122; pref=32*np.pi**2
def solveL(c): return brentq(lambda L: np.log(pref)-4*np.log(L)-c*L**2-np.log(Lam),1.0001,1e6)
for c,lab in [(3*np.pi**2/2,"round S^1xS^3"),(1.0,"O(1) instanton"),(0.2,"Higgs-critical")]:
    L=solveL(c)
    print(f"   c = {c:>6.2f} ({lab:>16}):  L = {L:>6.2f},  E = {1.22e19/L:.2e} GeV,  S = {c*L**2:.0f}")
