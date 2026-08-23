"""Scale hierarchy for the coarse-grained variation (Sec. IV.C of the paper).

Computes the defect separation l_def = n^{-1/4} from Lambda_obs, compares it
with the dark-energy length rho_Lambda^{-1/4}, and tabulates the Poisson
fluctuation dLambda/Lambda = (l_def/L)^2 versus coarse-graining scale.
"""
import numpy as np
lP=1.616e-35; alphabar=1.0; dchi=2.0; Lam=1e-122
n=Lam/(16*np.pi**2*alphabar*dchi)
l_def=(1/n)**0.25
rho=Lam/(8*np.pi); l_DE=(1/rho)**0.25
print(f"n lP^4          = {n:.3e}")
print(f"l_def           = {l_def:.3e} lP = {l_def*lP*1e3:.4f} mm")
print(f"rho_Lam^(-1/4)  = {l_DE:.3e} lP = {l_DE*lP*1e3:.4f} mm")
print(f"ratio           = {l_def/l_DE:.3f}\n")
print(f"{'L':>10} {'N = n V':>16} {'dLambda/Lambda':>18}")
for Lm,lab in [(1e-6,'1 um'),(1e-4,'0.1 mm'),(1e-3,'1 mm'),(1e-2,'1 cm'),
               (1.0,'1 m'),(1e3,'1 km'),(1.5e11,'1 AU')]:
    N=n*(Lm/lP)**4
    print(f"{lab:>10} {N:>16.3e} {1/np.sqrt(N):>18.3e}")
