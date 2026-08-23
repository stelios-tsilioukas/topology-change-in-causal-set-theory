import numpy as np
from flowtest import torus_pts, widest_betti

print("="*76)
print(" IS THE THRESHOLD SET BY HANDLE SIZE?  (vary the tube radius)")
print("="*76)
def torus_r(n,seed,R,r):
    rng=np.random.default_rng(seed)
    th=rng.uniform(0,2*np.pi,n); ph=rng.uniform(0,2*np.pi,n)
    return np.column_stack([(R+r*np.cos(th))*np.cos(ph),(R+r*np.cos(th))*np.sin(ph),r*np.sin(th)])

print(f"{'tube r':>8} {'N at which detection fails':>30}")
for r in (0.8,0.5,0.3):
    fail=None
    for N in (400,200,120,80,60,45,35,25,18):
        ok=0
        for s in range(10):
            P=torus_r(N,s,2.0,r)
            b=widest_betti(P)
            if b and b[1]>=1: ok+=1
        if ok/10 < 0.5: fail=N; break
    print(f"{r:>8.1f} {str(fail):>30}")

print("""
 The failure point tracks the POINT DENSITY needed to resolve the tube, i.e.
 the handle's physical size -- not any fixed combinatorial pattern size.

 CONCLUSION: a handle is not destroyed by removing 'its' n_h elements.  It is
 redundantly witnessed by very many sub-patterns, so the survival probability
 stays ~1 until the resolution exceeds the handle size, then collapses.
 The decimation law  p_h -> f^(n_h-1) p_h  is therefore NOT the behaviour of a
 physical topological feature.
""")
print("="*76)
print(" WHAT THE CORRECT FLOW LOOKS LIKE")
print("="*76)
print("""   n_resolved(l_cg) = density of handles of physical size  >~ l_cg
                    = integral over the foam size spectrum n(l_w)

   * all foam Planckian (l_w ~ l_P):  nothing is resolved for l_cg > l_P
                                      =>  dchi/dV|_cg = 0,  Lambda_eff = 0
   * scale-invariant foam:            n(>l_cg) ~ l_cg^-4  =>  Lambda ~ 1/l_cg^2
                                      matching Lambda_obs gives l_cg ~ 1/sqrt(Lambda)
                                      = the Hubble radius (everpresent Lambda,
                                      already excluded in Sec. VII)

   Either way the scale is set by the foam SIZE SPECTRUM, which causal set
   theory does not supply -- it is precisely the free input rho_w of the
   original continuum treatment.""")
