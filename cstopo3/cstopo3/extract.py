"""
cstopo3.extract -- order-only topology extraction in 3+1 dimensions.

Everything here touches the spacetime only through `st.precedes(i, js)`.  The
embedding is used nowhere below this line, which is the point of the exercise.

Pipeline
    antichain(st, t0, w)            inextendible antichain in a thin band
    shadows(st, A, cap_max)         past-shadows of the witnesses above it
    facets_at_cap(sh, cap)          Dowker facets at resolution `cap`
    scan(st, ...)                   Betti vector versus cap -- the plateau
"""
from __future__ import annotations
import numpy as np
import time
from .z2 import betti, faces_by_dim, collapse

__all__ = ["antichain", "shadows", "facets_at_cap", "scan", "witness_window"]


# ======================================================================
#  antichain
# ======================================================================
def antichain(st, t0=None, width=None, max_size=None, seed=0, verbose=False):
    """Greedy maximal antichain seeded from a thin time band.

    Two elements are unrelated iff |dt| <= d_spatial, so for a band of width
    comparable to the discreteness scale almost every pair in the band is
    unrelated and |A| approaches the band population, which scales as N^{3/4}.
    Candidates are considered in order of proximity to t0.
    """
    if t0 is None:
        t0 = 0.5 * st.T if not np.isfinite(st.t_split) else st.t_split - 3 * st.ell
    if width is None:
        width = 2.0 * st.ell
    cand = np.flatnonzero(np.abs(st.t - t0) <= 0.5 * width)
    if cand.size == 0:
        return np.array([], dtype=int)
    cand = cand[np.argsort(np.abs(st.t[cand] - t0))]
    aliveM = np.ones(cand.size, dtype=bool)
    A = []
    for k in range(cand.size):
        if not aliveM[k]:
            continue
        i = int(cand[k])
        A.append(i)
        rest = np.flatnonzero(aliveM)
        rest = rest[rest > k]
        if rest.size == 0:
            continue
        js = cand[rest]
        d = st.sdist(i, js)
        related = np.abs(st.t[js] - st.t[i]) > d
        aliveM[rest[related]] = False
        aliveM[k] = False
    A = np.array(sorted(A), dtype=int)
    if max_size is not None and A.size > max_size:
        rng = np.random.default_rng(seed)
        A = np.array(sorted(rng.choice(A, max_size, replace=False)))
    if verbose:
        print(f"    antichain: |A| = {A.size} from {cand.size} band candidates "
              f"(band {width:.4f}, ell {st.ell:.4f})")
    return A


# ======================================================================
#  witnesses and shadows
# ======================================================================
def witness_window(st, A, cap_max, margin=1.6):
    """Height above the antichain within which shadows can have <= cap_max
    elements.  |sigma(y)| ~ rho_A * (4/3)pi h^3, so h ~ (3 cap V / 4 pi |A|)^{1/3}."""
    vol = st.sa.volume
    h = (3.0 * cap_max * vol / (4.0 * np.pi * max(A.size, 1))) ** (1.0 / 3.0)
    return margin * h


def shadows(st, A, cap_max=16, tA=None, height=None, chunk=200_000,
            verbose=False):
    """Past-shadows sigma(y) = {a in A : a < y} for witnesses y above A.

    Returns (shadow_lists, n_witness, height) where shadow_lists is a list of
    sorted tuples of *positions in A* (not element indices).
    """
    if tA is None:
        tA = float(st.t[A].mean()) if A.size else 0.0
    if height is None:
        height = witness_window(st, A, cap_max)
    lo, hi = tA, tA + height
    W = np.flatnonzero((st.t > lo) & (st.t <= hi))
    if W.size == 0 or A.size == 0:
        return [], 0, height
    wpos = {int(y): k for k, y in enumerate(W)}
    acc = [[] for _ in range(W.size)]
    for k, a in enumerate(A):
        rel = st.precedes(int(a), W)
        for y in np.flatnonzero(rel):
            acc[y].append(k)
    out = [tuple(v) for v in acc if v]
    if verbose:
        sz = np.array([len(v) for v in out]) if out else np.array([0])
        print(f"    shadows: {W.size} witnesses in h={height:.4f}, "
              f"{len(out)} nonempty, |sigma| median {np.median(sz):.0f} "
              f"max {sz.max()}")
    return out, int(W.size), height


def facets_at_cap(shadow_lists, cap, nA, include_isolated=False):
    """Dowker facets at resolution `cap`, plus coverage diagnostics."""
    keep = set(s for s in shadow_lists if 0 < len(s) <= cap)
    covered = set()
    for s in keep:
        covered.update(s)
    facets = [s for s in keep]
    n_iso = nA - len(covered)
    if include_isolated:
        facets += [(k,) for k in range(nA) if k not in covered]
    return facets, len(covered), n_iso


# ======================================================================
#  the resolution scan
# ======================================================================
def scan(st, A=None, caps=(6, 8, 10, 12, 14, 16), maxdim=2,
         include_isolated=False, size_limit=4_000_000, verbose=True,
         shadow_data=None):
    """Betti vector versus the shadow-cardinality resolution.

    maxdim=2 gives (b0,b1,b2) and requires the complex to dimension 3.
    maxdim=1 gives (b0,b1) only and is roughly an order of magnitude cheaper;
    for a closed orientable 3-manifold b2 = b1 by Poincare duality, so b1 is
    already discriminating and b2 serves as a check on manifoldness.
    """
    if A is None:
        A = antichain(st, verbose=verbose)
    if shadow_data is None:
        shadow_data = shadows(st, A, cap_max=max(caps), verbose=verbose)
    sh, nW, h = shadow_data
    rows = []
    for c in caps:
        facets, ncov, niso = facets_at_cap(sh, c, A.size,
                                           include_isolated=include_isolated)
        if not facets:
            rows.append(dict(cap=c, betti=None, note="no facets",
                             n_facets=0, covered=ncov, isolated=niso))
            continue
        t0 = time.time()
        raw = faces_by_dim(facets, maxdim + 1)
        nraw = sum(len(v) for v in raw.values())
        if nraw > size_limit:
            rows.append(dict(cap=c, betti=None, note=f"too_big({nraw})",
                             n_facets=len(facets), covered=ncov,
                             isolated=niso, raw_total=nraw,
                             t=time.time() - t0))
            if verbose:
                print(f"    cap={c:3d}  SKIP  {nraw} raw simplices > limit")
            continue
        red = collapse(raw)
        b = betti(facets, maxdim=maxdim)
        dt = time.time() - t0
        rows.append(dict(cap=c, betti=b, note="", n_facets=len(facets),
                         covered=ncov, isolated=niso, raw_total=nraw,
                         collapsed_total=sum(len(v) for v in red.values()),
                         max_facet=max(len(f) for f in facets), t=dt))
        if verbose:
            print(f"    cap={c:3d}  beta={b}  facets={len(facets)} "
                  f"raw={nraw} coll={sum(len(v) for v in red.values())} "
                  f"cov={ncov}/{A.size} iso={niso}  {dt:.1f}s")
    return dict(A=A.size, n_witness=nW, height=h, ell=st.ell, rows=rows)


def plateau(rows):
    """The Betti vector taken over the largest number of consecutive caps."""
    vals = [tuple(r["betti"]) for r in rows if r["betti"] is not None]
    if not vals:
        return None
    from collections import Counter
    return Counter(vals).most_common(1)[0][0]

# ======================================================================
#  landmarks: the move that makes 3+1D tractable
# ======================================================================
def landmarks(A, nL, seed=0):
    """A sparse subset of the antichain, used as the vertex set of the complex.

    Why this is necessary.  For a point cloud sampling a 3-manifold the
    Dowker complex reproduces the topology only when the cover balls have
    radius r ~ 2 x the vertex spacing s, and the facet size is then
    m = (4pi/3)(r/s)^3 ~ 25-40 whatever the density.  Building that complex on
    the full antichain (|A| ~ 10^3-10^4) is hopeless; building it on a sparse
    subset L while keeping ALL of A as witnesses costs nothing in resolution,
    because the witnesses are what set the geometric scale.  The choice is a
    plain random subset, so it introduces no metric notion.
    """
    rng = np.random.default_rng(seed)
    return np.array(sorted(rng.choice(A, min(nL, A.size), replace=False)))


def landmark_shadows(st, L, cap, tA, volume, margin=1.35):
    """Shadows of the witnesses above tA, restricted to the landmark set.

    Returns (facets, n_witness, height, coverage).  A witness at height h has
    |sigma(y) & L| ~ rho_L (4pi/3) h^3, so the window that admits shadows up to
    `cap` is h ~ (3 cap V / 4 pi |L|)^{1/3}.
    """
    h = margin * (3.0 * cap * volume / (4.0 * np.pi * L.size)) ** (1.0 / 3.0)
    W = np.flatnonzero((st.t > tA) & (st.t <= tA + h))
    acc = [[] for _ in range(W.size)]
    for k, a in enumerate(L):
        for y in np.flatnonzero(st.precedes(int(a), W)):
            acc[y].append(k)
    facets = [tuple(v) for v in acc if 3 <= len(v) <= cap]
    cov = len({v for s in facets for v in s})
    return facets, int(W.size), h, cov


def slice_topology(st, t0, nL, cap, seed=0, volume=None, maxdim=1,
                   margin=1.35, verbose=False):
    """Betti numbers of the spatial slice at t0, from the causal order alone.

    maxdim=1 -> (b0,b1) by the cycle-space method: seconds.
    maxdim=2 -> (b0,b1,b2) by full Z2 reduction: minutes, and b2 needs a larger
                cap to converge (see run_b2_convergence.py).
    """
    from .fast01 import betti01
    from .z2 import betti as betti_full
    A = antichain(st, t0=t0)
    if A.size < nL:
        return None, dict(reason=f"|A|={A.size} < nL={nL}")
    L = landmarks(A, nL, seed=seed)
    if volume is None:
        volume = st.sa.volume if t0 >= st.t_split else st.sb.volume
    tA = float(st.t[A].mean())
    facets, nW, h, cov = landmark_shadows(st, L, cap, tA, volume, margin)
    if not facets:
        return None, dict(reason="no facets", n_witness=nW, height=h)
    info = dict(nA=int(A.size), nL=int(L.size), n_witness=nW, height=float(h),
                coverage=cov, n_facets=len(facets), ell=float(st.ell),
                max_facet=max(len(f) for f in facets))
    t1 = time.time()
    if maxdim == 1:
        b, extra = betti01(facets, L.size)
        info.update(extra)
    else:
        b = betti_full(facets, maxdim=maxdim)
    info["t_homology"] = time.time() - t1
    if verbose:
        print(f"      t0={t0:.2f} beta={b} |A|={A.size} |L|={L.size} "
              f"cov={cov}/{L.size} facets={len(facets)} h={h:.3f} "
              f"{info['t_homology']:.1f}s")
    return tuple(b), info
