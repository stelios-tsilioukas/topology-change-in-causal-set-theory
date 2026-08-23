import math
pi=math.pi
lP2=1.0
Lam_obs=2.85e-122     # Lambda lP^2, Planck 2018
S1=287.0              # paper's wormhole action

print("="*72); print("1. NARIAI ACTION, IF ONE GRANTS IT")
print("="*72)
S_nar=-2*pi/Lam_obs
print(f"  S_Nariai = -2pi/(G Lambda) = {S_nar:.3e}   (NEGATIVE)")
print(f"  weight e^-S = e^(+{-S_nar:.2e})  -- conformal-factor problem of Sec V.E, extreme form")
print(f"  compare S^4 (de Sitter): {-3*pi/Lam_obs:.3e}")
print("  and it scales as 1/Lambda: the species is defined only if Lambda is already there.")

print()
print("="*72); print("2. TOPOLOGICAL CLASSIFICATION OF DEFECTS")
print("="*72)
print("  closed orientable connected X:  chi = 2 - 2b1 + b2   (b3=b1, b4=b0=1)")
print("  glue in as a handle (remove two 4-balls):  d_chi = chi(X) - 2 = b2 - 2*b1")
print("  axion charge  n = int_{S^3} H_3  requires  b3 = b1 >= 1")
print()
rows=[("S^1 x S^3",1,0),("S^2 x S^2",0,2),("S^4",0,0),("CP^2",0,1),("K3",0,22),
      ("(S^1xS^3)#(S^2xS^2)",1,2),("(S^1xS^3)#2(S^2xS^2)",1,4),("T^4",4,6)]
print(f"  {'X':<22}{'b1':>4}{'b2':>4}{'chi':>6}{'d_chi':>7}   axion charge?")
for nm,b1,b2 in rows:
    chi=2-2*b1+b2
    print(f"  {nm:<22}{b1:>4}{b2:>4}{chi:>6}{chi-2:>7}   {'YES' if b1>=1 else 'no'}")
print()
print("  => simply connected  =>  b1=b3=0  =>  NO axion charge.  Kills S^2xS^2, S^4, CP^2, K3.")
print("  => minimal charged defect: b1=1,b2=0 = S^1 x S^3, d_chi = -2  <- the GS wormhole")
print("  => a CHARGED defect with d_chi>=0 needs b2 >= 2b1 >= 2, i.e. it must contain")
print("     the very S^2xS^2 factors shown in part 1 to require Lambda>0 to be supported.")

print()
print("="*72); print("3. HOW MUCH ROOM IS LEFT? (falsifiable bound)")
print("="*72)
print("  Lambda_eff > 0  <=>  sum_i d_chi_i n_i < 0  <=>  2A_-e^-S_- > 2A_+e^-S_+")
print(f"  with comparable prefactors:   S_+ > S_- = {S1:.0f}")
print("  i.e. the sign is safe unless some d_chi>0 saddle has Euclidean action BELOW 287.")
for dS in (0,5,10,25,50):
    print(f"     S_+ = {S1+dS:5.0f}:  Nariai-sector contamination = e^-{dS} = {math.exp(-dS):.2e}")

print()
print("="*72); print("4. THE GB TERM SELF-WEIGHTS THE FOAM (new, secondary)")
print("="*72)
print("  Each defect's own action contains the GB piece:")
print("    |S_GB| = 16 pi^2 alpha |d_chi| / kappa^2 = 2 pi abar |d_chi|")
print("  so the two species are split by  exp(8 pi abar)  BEFORE any gravitational difference:")
for ab in (1e-4,1e-2,1,1e2,1e4):
    print(f"     abar={ab:8.0e}:  2pi*abar*|d_chi|={4*pi*ab:12.4g}   split e^(8 pi abar) = {'%.2e'%math.exp(min(8*pi*ab,700))}")
print()
print("  Consequence for Sec V.G: with the self-weight included the matching becomes")
print("     S_grav = 287.0 + ln(abar) + 4*pi*abar     (sign of last term = Euclidean convention)")
print(f"  {'abar':>10}{'S_grav':>12}{'f [GeV]':>14}{'vs paper':>12}")
for ab in (1e-4,1e-2,1,1e2,1e4):
    Sg=287.0+math.log(ab)+4*pi*ab
    f=(math.sqrt(6)*pi/2)*2.435e18/Sg
    f0=(math.sqrt(6)*pi/2)*2.435e18/(287.0+math.log(ab))
    print(f"  {ab:>10.0e}{Sg:>12.1f}{f:>14.3e}{f/f0:>11.3f}x")
print("  => the claimed 6%-over-eight-decades stability holds only for abar <~ 0.1.")
