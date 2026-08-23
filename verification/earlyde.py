import numpy as np

print("="*70)
print(" EARLY DARK ENERGY TEST FOR H^2-TRACKING TOPOLOGICAL DARK ENERGY")
print("="*70)

# ---------------------------------------------------------------
# Step 0: does Omega_DE really stay constant?  Check V ~ H^-4 per era
# ---------------------------------------------------------------
print("\n[0] Is V ~ H^-4 robust across eras?  (a ~ t^n)")
print(f"{'era':>12} {'n':>6} {'V(t) scaling':>16} {'H^-4 scaling':>16}  same?")
for era,n in [("radiation",0.5),("matter",2.0/3.0)]:
    # horizon d_H = t/(1-n); V ~ int d_H^3 dt ~ t^4/(4(1-n)^3);  H = n/t -> H^-4 = t^4/n^4
    cV = 1.0/(4*(1-n)**3); cH = 1.0/n**4
    print(f"{era:>12} {n:>6.3f} {'t^4 x %.2f'%cV:>16} {'t^4 x %.2f'%cH:>16}   YES (both t^4)")
print("  => V ∝ H^-4 in any single-component era, so sigma_Lambda ∝ H^2 exactly,")
print("     hence Omega_DE = sigma_Lambda/(3H^2) = CONSTANT at all epochs.")

# ---------------------------------------------------------------
# Step 1: fix sqrt(p) from today's dark energy
# ---------------------------------------------------------------
print("\n[1] Calibrate p from Omega_DE,0 = 0.70")
Om_DE0 = 0.70
# Omega_DE = Lambda/(3H^2);  Lambda = 32 pi^2 sqrt(p) H^2  =>  Omega = 32 pi^2 sqrt(p)/3
sqrtp = 3*Om_DE0/(32*np.pi**2)
p = sqrtp**2
print(f"    Omega_DE = 32 pi^2 sqrt(p)/3   =>   sqrt(p) = {sqrtp:.4e},  p = {p:.3e}")
print(f"    (combinatorial estimate from Sec.IV was p ~ 1e-4, range 1e-1..1e-7)")
print(f"    -> CONSISTENT: required p sits inside the predicted range.")

# ---------------------------------------------------------------
# Step 2: what does constant Omega_DE do at BBN?
# ---------------------------------------------------------------
print("\n[2] Consequence at BBN:  H is boosted by 1/sqrt(1-Omega_DE)")
def boost(x): return 1.0/np.sqrt(1.0-x)
def Cr(x):    return (boost(x)-1.0)**2        # the C_r metric of arXiv:2507.18389
def dNeff(x):
    # rho_extra/rho_rad = x/(1-x);  rho_rad = rho_gamma[1+0.2271 Neff]
    f = 1.0 + 0.2271*3.046
    return f*(x/(1.0-x))/0.2271

print(f"    Omega_DE = {Om_DE0}:  H/H_std = {boost(Om_DE0):.3f}  ({100*(boost(Om_DE0)-1):.0f}% faster)")
print(f"                          C_r     = {Cr(Om_DE0):.3f}")
print(f"                          dN_eff  = {dNeff(Om_DE0):.1f}")

# ---------------------------------------------------------------
# Step 3: bounds
# ---------------------------------------------------------------
print("\n[3] Observational bounds vs prediction")
bounds = [("BBN  (dN_eff < 0.4)",  "dNeff", 0.4),
          ("BBN  (C_r  < 0.01)",   "Cr",    0.01),
          ("CMB  (Omega_EDE<0.02)","Om",    0.02)]
from scipy.optimize import brentq
print(f"{'bound':>24} {'max Omega_DE':>14} {'max sqrt(p)':>13} {'exceeded by':>13}")
rows=[]
for name,kind,val in bounds:
    if kind=="dNeff":  xmax = brentq(lambda x: dNeff(x)-val, 1e-8, 0.999)
    elif kind=="Cr":   xmax = brentq(lambda x: Cr(x)-val,    1e-8, 0.999)
    else:              xmax = val
    sp = 3*xmax/(32*np.pi**2)
    rows.append((name,xmax,sp))
    print(f"{name:>24} {xmax:>14.4f} {sp:>13.3e} {Om_DE0/xmax:>12.1f}x")

# ---------------------------------------------------------------
# Step 4: the irreducible tension
# ---------------------------------------------------------------
print("\n[4] The tension")
xmax_all = min(r[1] for r in rows)
print(f"    Tightest bound:      Omega_DE < {xmax_all:.4f}")
print(f"    Required today:      Omega_DE = {Om_DE0}")
print(f"    Ratio:               {Om_DE0/xmax_all:.0f}x  (in p: {(Om_DE0/xmax_all)**2:.0f}x)")
print(f"    If we SATISFY the bound, today's DE fraction is only {xmax_all:.3f},")
print(f"    i.e. the model supplies {100*xmax_all/Om_DE0:.1f}% of the observed dark energy.")

# ---------------------------------------------------------------
# Step 5: how much sign-flip averaging could rescue it?
# ---------------------------------------------------------------
print("\n[5] Does mean-zero sign flipping rescue BBN?")
print("    Lambda flips sign ~once per Hubble time; BBN spans t ~ 1s -> 1000s,")
print("    i.e. ~ln(1000) = %.1f e-folds  =>  ~%d independent Hubble patches in time."%(np.log(1000),int(np.log(1000))))
N_flip = np.log(1000)
print(f"    Random-walk suppression of the INTEGRATED effect: ~1/sqrt(N) = {1/np.sqrt(N_flip):.2f}")
print(f"    => effective |Omega_DE| during BBN ~ {Om_DE0/np.sqrt(N_flip):.2f}, still >> {xmax_all:.3f}")
print("    Sign flipping does NOT rescue it: the INSTANTANEOUS |Omega_DE| ~ 0.7")
print("    modulates H by +/-40% throughout, which Y_p cannot tolerate.")

# ---------------------------------------------------------------
# Step 6: what suppression is needed
# ---------------------------------------------------------------
print("\n[6] Required epoch-dependence of Gamma (open problem (i))")
print("    To satisfy BBN and still give Omega_DE,0 = 0.70, need Gamma_i(z) with")
z_bbn = 1e9; z_rec=1100
for nm,z,xm in [("BBN", z_bbn, xmax_all),("recomb", z_rec, 0.02)]:
    supp = (xm/Om_DE0)
    # sigma ~ sqrt(p) ; Omega ∝ sqrt(p)  => p must be suppressed by supp^2
    n_index = np.log(1/supp)/np.log(1+z)
    print(f"      {nm:>7}: suppression of sqrt(p) by {supp:.4f}  (p by {supp**2:.1e})")
    print(f"               if p(z) ∝ (1+z)^-m  =>  m = {2*n_index:.3f}")
