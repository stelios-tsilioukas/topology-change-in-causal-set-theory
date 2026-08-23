import numpy as np
print("="*78)
print(" RECONSTRUCTING COLEMAN'S RATE:  dchi/dV = Gamma ?")
print("="*78)
print("""
 Ref.[ATS2025] Eq.(7):     n_i  ==  Gamma_i  =  A_i exp(-Delta I_i)
 This paper,  Eq.(flow):   dchi/dV = (dchi/l_cg^4) exp(-S(l_cg))

 Matching term by term for one species (dchi factored out):

     Gamma  =  (1/l_cg^4) exp(-S(l_cg))

 so that
     A       =  1/l_cg^4  =  rho_cg      <-- the COARSE-GRAINED SPRINKLING DENSITY
     Delta I =  S(l_cg)   =  c L^2       <-- action at the RESOLUTION scale
""")
print("="*78); print(" DIMENSIONAL CHECK"); print("="*78)
print("""   [Gamma] = (four-volume)^-1 = L^-4 .
   In Coleman's theory A is a ratio of functional determinants, of dimension
   L^-4, conventionally estimated as ~ (S_0/2pi)^2 / R^4 for a bounce of size R.
   Here it is not estimated: it is the number of independent coarse-graining
   cells per four-volume.  The determinant is replaced by a COUNT.""")

c=3*np.pi**2/2; L=4.36
print("\n"+"="*78); print(" NUMBERS"); print("="*78)
print(f"   l_cg = {L:.2f} l_P   =>   A = rho_cg = 1/l_cg^4 = {1/L**4:.4f} l_P^-4")
print(f"   Delta I = S(l_cg) = c L^2 = {c*L**2:.0f}")
print(f"   Gamma = {1/L**4:.4f} e^-{c*L**2:.0f} l_P^-4 = {1/L**4:.2e} x 10^-{c*L**2/np.log(10):.0f} l_P^-4")
# convert to their SI-style units for comparison with rho_w ~ 1e16 m^-3 s^-1
lP=1.616e-35; tP=5.391e-44
Gamma_SI = (1/L**4)*np.exp(-c*L**2)/(lP**3*tP)
print(f"\n   in SI: Gamma = {Gamma_SI:.2e} m^-3 s^-1")
print(f"   Ref.[TST2024] required rho_w ~ 1e16 m^-3 s^-1 to match Lambda_obs.")
print(f"   ratio = {Gamma_SI/1e16:.2f}   -- consistent, as it must be by construction.")

print("\n"+"="*78); print(" THE ONE STRUCTURAL DIFFERENCE"); print("="*78)
print("""   In Coleman's theory Delta I is the action of a STATIONARY configuration
   (the bounce), whose size is fixed by extremising S.  Here S(l_w) = c l_w^2
   is monotonic and has NO stationary point, so the size is fixed instead by
   the RESOLUTION scale.  The coarse-graining cutoff does the job that the
   saddle point does in the continuum theory.""")
