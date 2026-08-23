import math
from scipy.optimize import brentq
pi=math.pi; lP_mm=1.616255e-32; Mp=2.435e18

# --- observed Lambda, Planck 2018 ---
H0=67.4*1000/3.0857e22; c_=2.99792458e8; lP=1.616255e-35
Lam=3*0.685*(H0/c_)**2; Lobs=Lam*lP**2
print(f"Lambda = {Lam:.4e} m^-2 ;  Lambda lP^2 = {Lobs:.4e}")

cc=4/(3*pi**2); K=8/cc**2
print(f"a0^2/lP^2 = {cc:.5f} S ;  K = {K:.4f}")

print("\n--- MAIN TABLE (prefactor (S/2pi)^2 a0^-4) ---")
for ab in (1e-4,1e-2,1,1e2,1e4):
    S=math.log(K*ab/Lobs); f=(math.sqrt(6)*pi/2)*Mp/S
    print(f"  abar={ab:8.0e}  S1={S:7.2f}  f={f:.4e} GeV  a0/lP={math.sqrt(cc*S):.3f}")
S1=math.log(K/Lobs); print(f"  => abar=1: S1={S1:.2f}, f={(math.sqrt(6)*pi/2)*Mp/S1:.4e}, a0={math.sqrt(cc*S1):.3f}")

print("\n--- ALTERNATIVE prefactor (no S^2, i.e. eq:spectrum's own) ---")
g=lambda S: math.log(32*pi**2/(cc*S)**2)-S-math.log(Lobs)
Sa=brentq(g,100,600); fa=(math.sqrt(6)*pi/2)*Mp/Sa
f1=(math.sqrt(6)*pi/2)*Mp/S1
print(f"  S1={Sa:.2f}  f={fa:.4e} GeV   spread vs main = {abs(fa-f1)/f1*100:.1f}%")

print("\n--- defect density and scales ---")
nbar=Lobs/(32*pi**2); ldef=nbar**-0.25
print(f"  n lP^4 = {nbar:.4e}   l_def = {ldef:.4e} lP = {ldef*lP_mm:.4f} mm")
rho=Lobs/(8*pi); lrho=rho**-0.25
print(f"  rho=Lam/8piG: rho^-1/4 = {lrho*lP_mm:.4f} mm")
hbarc=1.97327e-7
rho_phys=Lam*c_**2/(8*pi*6.674e-11)
E=(rho_phys*c_**2*(1.0545718e-34*c_)**3)**0.25/1.602e-19
print(f"  physical rho_Lambda^1/4 = {E*1e3:.3f} meV -> {hbarc/E*1e3:.4f} mm  (agrees)")
print(f"  RATIO l_def/rho^-1/4 = {ldef/lrho:.4f}   and (4pi)^(1/4) = {(4*pi)**0.25:.4f}  <-- EXACT identity")

print("\n--- fluctuation table ---")
for L,lab in [(0.1,'0.1 mm'),(1,'1 mm'),(10,'1 cm'),(1000,'1 m')]:
    N=(L/(ldef*lP_mm))**4
    print(f"  L={lab:6s}  N={N:.3g}   dL/L={1/math.sqrt(N):.3g}")

print("\n--- round-handle resolvability check (v5 route), c=3pi^2/2 ---")
cR=3*pi**2/2
h=lambda L: math.log(32*pi**2)-4*math.log(L)-cR*L*L-math.log(Lobs)
Lr=brentq(h,1,50); print(f"  c={cR:.3f}  l_cg={Lr:.3f} lP   S(l_cg)={cR*Lr*Lr:.1f}")
for ab in (1e-2,1e2):
    hh=lambda L: math.log(32*pi**2*ab)-4*math.log(L)-cR*L*L-math.log(Lobs)
    print(f"    abar={ab:.0e}: l_cg={brentq(hh,1,50):.3f} lP")

print("\n--- GB self-weight ---")
print(f"  |S_GB| per defect = 2 pi abar |dchi| = {4*pi:.3f} abar ; split e^(8 pi abar)")
for ab in (1e-2,0.1,1):
    print(f"    abar={ab:5.2f}:  4 pi abar = {4*pi*ab:8.3f}  -> shifts S1 by this much")
print("\n--- renormalisation violation ---")
print(f"  e^-S1 = Lobs/(K abar) = {Lobs/K:.2e}/abar")
