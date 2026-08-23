import numpy as np
print("="*72); print(" COARSE-GRAINING FLOW OF THE TOPOLOGICAL DENSITY"); print("="*72)
print("""
 Decimation: keep each element with probability f.  Sprinkling -> sprinkling
 at density f*rho_c, so the effective discreteness scale is
      l_cg = l * f^(-1/4)        <=>   f = (l/l_cg)^4 = L^-4 ,  L = l_cg/l

 A handle occupies n_h elements.  Survival requires all n_h to survive:
      N_h -> f^{n_h} N_h ,   N -> f N
      p_h -> f^{n_h - 1} p_h                      (density per element)

 The four-volume density then flows as
      dchi/dV|_cg = dchi * p_h(f) / l_cg^4
                  = (dchi * p_h / l^4) * f^{n_h-1} * f
                  = (dchi * p_h / l^4) * L^{-4 n_h}
""")
nh=6
print(f"   with n_h = {nh}:   dchi/dV  ∝  L^(-{4*nh})\n")

print("="*72); print(" WHAT SCALE REPRODUCES THE OBSERVED LAMBDA?"); print("="*72)
Lam_obs=1e-122          # Lambda * l_P^2
pref=32*np.pi**2        # 16 pi^2 * |dchi| with |dchi|=2, alpha = l_P^2
print(f"   Lambda*l^2 = 32 pi^2 * p_h * L^(-{4*nh})  =  {Lam_obs:.0e}\n")
print(f"{'p_h':>10} {'L = l_cg/l_P':>16} {'l_cg [m]':>13} {'E = M_Pl/L [GeV]':>18}")
for p in [1e-1,1e-2,1e-4,1e-6,1e-7]:
    L=(pref*p/Lam_obs)**(1.0/(4*nh))
    lcg=1.616e-35*L
    E=1.22e19/L
    print(f"{p:>10.0e} {L:>16.3e} {lcg:>13.2e} {E:>18.2e}")
print("""
   Remarkably insensitive: p_h spanning 6 orders of magnitude moves L by <2x,
   because it enters as p^(1/24).  The scale is essentially FIXED:
        l_cg ~ 10^5 l_P  ~  10^-30 m  ~  10^14 GeV
   i.e. the GUT / seesaw / upper-inflation-scale window.
""")
print("="*72); print(" CONSEQUENCE FOR w"); print("="*72)
print("""   Manifoldlikeness makes the LOCAL abundances universal, hence p_h(l_cg)
   epoch-INDEPENDENT.  In the MEAN channel that gives

        Lambda_eff = -16 pi^2 alpha <dchi/dV> = CONSTANT   =>   w = -1 exactly.

   The same universality that killed the fluctuation channel (b=0 => w=0)
   *guarantees* w = -1 in the mean channel.""")
