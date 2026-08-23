import numpy as np
from scipy.optimize import brentq

def percolate(n, p, rng):
    """Transitive percolation: relate i<j w.p. p, then transitive closure."""
    R = np.triu(rng.random((n,n)) < p, 1)
    # transitive closure (Warshall, vectorised)
    for k in range(n):
        R |= np.outer(R[:,k], R[k,:])
    return np.triu(R,1)

def ordering_fraction(R):
    n = R.shape[0]
    return R.sum()/(n*(n-1)/2)

def interval_abundances(R, kmax=4):
    """n_k = (# k-element order intervals)/N  -- the Benincasa-Dowker observables.
       interval(i,j) = {m : i<m<j}; count intervals of cardinality k-2 interior."""
    n = R.shape[0]
    counts = np.zeros(kmax+1)
    idx = np.argwhere(R)
    for i,j in idx:
        card = int((R[i,:] & R[:,j]).sum())   # interior size
        if card+2 <= kmax: counts[card+2] += 1
    return counts/n

def mean_r(n, p, seeds=6):
    return np.mean([ordering_fraction(percolate(n,p,np.random.default_rng(s))) for s in range(seeds)])

print("="*72)
print(" CSG COSMIC RENORMALISATION:  how must p run as the causet grows?")
print("="*72)
print("\nCriterion: a manifoldlike causet sprinkled into a fixed-shape region has a")
print("FIXED ordering fraction r (a purely geometric number). Demand r = const as N grows.\n")

r_target = 0.30
Ns = [30, 50, 80, 120, 180, 260]
print(f"{'N':>6} {'p*(N) for r=0.30':>20} {'r achieved':>12}")
ps=[]
for n in Ns:
    f = lambda pp: mean_r(n, pp) - r_target
    lo, hi = 1e-4, 0.9
    try:
        pstar = brentq(f, lo, hi, xtol=1e-5)
    except ValueError:
        pstar = np.nan
    ps.append(pstar)
    print(f"{n:>6} {pstar:>20.5f} {mean_r(n,pstar):>12.4f}")

ps=np.array(ps); Ns_a=np.array(Ns,float)
good = ~np.isnan(ps)
a, b = np.polyfit(np.log(Ns_a[good]), np.log(ps[good]), 1)
print(f"\n  FIT:  p*(N) ∝ N^({a:.3f})      [exponent a = {a:.3f}]")
print(f"  => p DECREASES with growth" if a<0 else "  => p INCREASES with growth")
