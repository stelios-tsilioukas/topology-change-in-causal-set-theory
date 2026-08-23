import numpy as np
from math import comb

# ---------- Sorkin's 4D smearing function ----------
def f4(n, eps):
    r = eps/(1.0-eps)
    return (1-eps)**n * (1 - 9*r*n + 8*r**2*n*(n-1) - (4.0/3.0)*r**3*n*(n-1)*(n-2))

# ---------- the layer coefficients of the sharp 4D BD operator ----------
C4 = [1.0, -9.0, 16.0, -8.0]     # c_k for k = 0,1,2,3 elements strictly between

# ---------- claim: sum_k c_k C(n,k) eps^k (1-eps)^{n-k}  ==  f4(n,eps) ----------
print("="*74)
print(" IDENTITY CHECK:  binomial thinning of the sharp layers == Sorkin smearing")
print("="*74)
print(f"{'n':>4} {'eps':>7} {'binomial thinning':>22} {'f4(n,eps)':>18} {'match':>8}")
ok=True
for eps in (0.05,0.2,0.5,0.8):
    for n in (0,1,2,3,5,9,17):
        thin = sum(C4[k]*comb(n,k)*eps**k*(1-eps)**(n-k) for k in range(4) if k<=n)
        val  = f4(n,eps)
        m = abs(thin-val) < 1e-10*max(1,abs(val))
        ok &= m
        if n in (0,2,9,17):
            print(f"{n:>4} {eps:>7.2f} {thin:>22.12f} {val:>18.12f} {str(m):>8}")
print(f"\n  ALL {'MATCH' if ok else 'FAIL'}  -> the identity is exact.\n")

print("  Analytic proof: decimation with retention eps sends an interval of")
print("  cardinality n to one of cardinality k with probability C(n,k) eps^k (1-eps)^(n-k).")
print("  Hence  <sum_k c_k>  =  sum_k c_k C(n,k) eps^k (1-eps)^(n-k)")
print("                      =  (1-eps)^n [1 - 9nr + 16 C(n,2) r^2 - 8 C(n,3) r^3],  r=eps/(1-eps)")
print("  and  16*C(n,2) = 8n(n-1),  8*C(n,3) = (4/3)n(n-1)(n-2)   ==>  exactly f4(n,eps).  QED")
