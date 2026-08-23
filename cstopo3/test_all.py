#!/usr/bin/env python3
"""
test_all.py -- run this first on new hardware.  Should take under two minutes.

Checks, in order of what they would break:
  1. the homology engines against exact triangulations with known Betti numbers,
     including the wedge S^1 v S^2 v S^3 which shares a Betti vector with
     S^1 x S^2 and is the negative control for the manifold condition;
  2. agreement between the fast (b0,b1) path and the general one;
  3. the geometry: volumes, exact-distance symmetry, sprinkling density
     continuity across t_split;
  4. |A| ~ N^{3/4};
  5. end-to-end recovery of beta(S^3) and beta(S^1 x S^2) from the order alone.
"""
from __future__ import annotations
import sys, os, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cstopo3.validate import KNOWN
from cstopo3.z2 import betti
from cstopo3.fast01 import betti01
from cstopo3.geom3 import S3, S1xS2, S3Warped, Surgery, Spacetime
from cstopo3.extract import antichain, slice_topology

FAIL = []


def check(name, cond, detail=""):
    print(f"  {'PASS' if cond else 'FAIL'}  {name}{('  ' + detail) if detail else ''}")
    if not cond:
        FAIL.append(name)


def relabel(facets):
    vs = sorted({v for f in facets for v in f})
    rl = {v: i for i, v in enumerate(vs)}
    return [tuple(sorted(rl[v] for v in f)) for f in facets], len(vs)


print("== 1. homology engine vs exact triangulations ==")
for name, (facets, expect) in KNOWN.items():
    b = betti(facets, maxdim=len(expect) - 1)
    check(f"betti {name}", b == expect, f"{b}")

print("\n== 2. fast (b0,b1) agrees with the general path ==")
for name, (facets, expect) in KNOWN.items():
    fx, nv = relabel(facets)
    (b0, b1), _ = betti01(fx, nv)
    check(f"betti01 {name}", (b0, b1) == expect[:2], f"({b0},{b1})")

print("\n== 3. geometry ==")
check("vol(S3) = 2 pi^2", abs(S3(1.0).volume - 2 * np.pi ** 2) < 1e-9)
check("vol(S1xS2) = 8 pi^2", abs(S1xS2(1.0, 1.0).volume - 8 * np.pi ** 2) < 1e-9)
check("vol(S3warped) = vol(S3)",
      abs(S3Warped(1.0, neck=1.5).volume - 2 * np.pi ** 2) < 1e-3,
      f"{S3Warped(1.0, neck=1.5).volume:.4f}")
sg = Surgery(1.0, 0.4, neck=1.5)
check("surgery 1-cycle length = 2*neck+pi",
      abs(sg.circumference - (2 * 1.5 + np.pi)) < 1e-9, f"{sg.circumference:.3f}")

rng = np.random.default_rng(0)
for sl in (S3(1.0), S1xS2(1.0, 1.0)):
    st = Spacetime(sl, 400, T=2.0, seed=0)
    d_ij = sl.dist(st.x, 5, np.array([9]))[0]
    d_ji = sl.dist(st.x, 9, np.array([5]))[0]
    check(f"{sl.name} distance symmetric", abs(d_ij - d_ji) < 1e-12)
    check(f"{sl.name} d(i,i) = 0", abs(sl.dist(st.x, 5, np.array([5]))[0]) < 1e-9)

st = Spacetime(S3Warped(1.0, neck=1.5), 20000, T=4.5, seed=1,
               slice_after=Surgery(1.0, 0.4, neck=1.5), t_split=2.0)
rho_b = st.pre.sum() / (S3Warped(1.0, neck=1.5).volume * 2.0)
rho_a = (~st.pre).sum() / (Surgery(1.0, 0.4, neck=1.5).volume * 2.5)
check("sprinkling density continuous across t_split",
      abs(rho_b - rho_a) / rho_b < 0.02, f"{rho_b:.1f} vs {rho_a:.1f}")

print("\n== 4. |A| scales as N^{3/4} ==")
sizes = []
for N in (5000, 20000, 80000):
    A = antichain(Spacetime(S3(1.0), N, T=3.0, seed=1), t0=1.2)
    sizes.append(A.size)
    print(f"       N={N:6d}  |A|={A.size:5d}")
p = np.polyfit(np.log([5000, 20000, 80000]), np.log(sizes), 1)[0]
check("exponent close to 0.75", abs(p - 0.75) < 0.12, f"fitted {p:.3f}")

print("\n== 5. end to end: beta from the order alone ==")
for sl, want in ((S3(1.0), (1, 0)), (S1xS2(1.0, 1.0), (1, 1))):
    t0 = time.time()
    st = Spacetime(sl, 20000, T=3.0, seed=1)
    b, info = slice_topology(st, 1.2, 250, 14, seed=1, volume=sl.volume, maxdim=1)
    check(f"{sl.name} -> {want}", b == want,
          f"got {b}, |A|={info.get('nA')}, {time.time()-t0:.1f}s")

print("\n" + ("ALL TESTS PASSED" if not FAIL else f"*** {len(FAIL)} FAILED: {FAIL}"))
sys.exit(1 if FAIL else 0)
