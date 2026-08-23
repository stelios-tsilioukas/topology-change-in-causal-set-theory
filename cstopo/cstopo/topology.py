"""
cstopo.topology -- order-theoretic topology extraction.

Everything here uses ONLY the causal matrix.  No embedding coordinates are
touched, which is the point: the demonstration must show that topology is
recoverable from the order alone.

Pipeline:
    antichain(R, ...)            inextendible antichain (a "moment of time")
    thicken(R, A, n)             elements within n causal layers to the future
    mrs_nerve(R, A, n)           maximal simplices of the past-shadow cover
    betti_gf2(maximal)           Betti numbers over GF(2)
    stability_scan(R, A, ns)     Betti vector vs thickening -- the plateau
"""
from __future__ import annotations
import numpy as np
from itertools import combinations

__all__ = ["antichain_from_band", "layer_index", "mrs_nerve", "betti_gf2",
           "euler_char", "stability_scan", "morse_count"]


# ======================================================================
#  antichains
# ======================================================================
def antichain_from_band(R, t, t0, width, max_size=None, rng=None):
    """Greedy inextendible antichain seeded from the time band |t - t0| < width/2.

    Elements are added in order of proximity to t0, skipping any that are
    causally related to one already chosen.  This yields an antichain that is
    inextendible within the band; for a thin band in a manifoldlike causet it
    approximates a spatial hypersurface.
    """
    cand = np.flatnonzero(np.abs(t - t0) <= width / 2)
    cand = cand[np.argsort(np.abs(t[cand] - t0))]
    A = []
    masks = []
    for i in cand:
        mi = R.row_mask(i)
        rel = False
        for a, ma in zip(A, masks):
            if mi[a] or ma[i]:
                rel = True
                break
        if not rel:
            A.append(int(i))
            masks.append(mi)
    A = np.array(sorted(A))
    if max_size is not None and len(A) > max_size:
        rng = rng or np.random.default_rng(0)
        A = np.array(sorted(rng.choice(A, max_size, replace=False)))
    return A


def layer_index(R, A, nmax=8):
    """Causal layer of each element above the antichain A.

    layer = 0 for elements of A; layer = k for elements whose longest chain
    from A has k links.  Elements not to the future of A get -1.
    """
    N = R.n
    layer = np.full(N, -1, dtype=np.int16)
    layer[A] = 0
    frontier = set(int(a) for a in A)
    for k in range(1, nmax + 1):
        nxt = set()
        for i in frontier:
            for j in R.row(i):
                j = int(j)
                if layer[j] < 0:
                    layer[j] = k
                    nxt.add(j)
        if not nxt:
            break
        frontier = nxt
    return layer


# ======================================================================
#  MRS nerve
# ======================================================================
def mrs_nerve(R, A, layer, nthick):
    """Maximal simplices of the nerve of the past-shadow cover.

    For each element y within nthick layers to the future of A, its shadow is
        sigma(y) = { a in A : a precedes y }.
    The nerve has A as vertices and the sigma(y) as (candidate) simplices; we
    return the maximal ones, plus isolated vertices.
    """
    slab = np.flatnonzero((layer >= 1) & (layer <= nthick))
    pos = {int(a): k for k, a in enumerate(A)}
    shadows = set()
    for y in slab:
        m = R.row_mask(int(y)) if False else None
        # collect a in A with a -> y
        s = tuple(sorted(pos[int(a)] for a in A if R.row_mask(int(a))[int(y)]))
        if s:
            shadows.add(s)
    mx = [s for s in shadows if not any(set(s) < set(u) for u in shadows)]
    covered = set()
    for s in mx:
        covered.update(s)
    for k in range(len(A)):
        if k not in covered:
            mx.append((k,))
    return mx


def mrs_nerve_fast(R, A, layer, nthick, shadow_cap=None):
    """Nerve of the past-shadow cover, thickening by shadow SIZE.

    The MRS prescription requires the cover sets to be *small*: an element far
    to the future of A has essentially all of A in its past, and its shadow
    carries no local information.  We therefore keep only witnesses whose
    shadow has cardinality <= shadow_cap, which is the causal-interval measure
    of "close to A" and is the correct discrete analogue of a small ball.

    nthick still restricts to the first nthick causal layers; shadow_cap is the
    resolution parameter whose plateau we scan.
    """
    slab = np.flatnonzero((layer >= 1) & (layer <= nthick))
    nA = len(A)
    if slab.size == 0:
        return [(k,) for k in range(nA)]
    slabset = np.zeros(R.n, dtype=bool)
    slabset[slab] = True
    acc = {}
    for k, a in enumerate(A):
        fut = R.row(int(a))
        fut = fut[slabset[fut]]
        for y in fut:
            acc.setdefault(int(y), []).append(k)
    if shadow_cap is None:
        shadow_cap = nA
    shadows = set(tuple(sorted(v)) for v in acc.values()
                  if 0 < len(v) <= shadow_cap)
    mx = [s for s in shadows if not any(set(s) < set(u) for u in shadows)]
    covered = set()
    for s in mx:
        covered.update(s)
    for k in range(nA):
        if k not in covered:
            mx.append((k,))
    return mx


# ======================================================================
#  homology over GF(2)  (no external dependencies)
# ======================================================================
def _all_faces(maximal, maxdim):
    faces = {}
    for s in maximal:
        s = tuple(sorted(s))
        top = min(len(s), maxdim + 1)
        for k in range(1, top + 1):
            for c in combinations(s, k):
                faces.setdefault(k - 1, set()).add(c)
    return {d: sorted(v) for d, v in faces.items()}


def _gf2_rank(cols):
    piv = {}
    rank = 0
    for v in cols:
        cur = v
        while cur:
            h = cur.bit_length() - 1
            if h in piv:
                cur ^= piv[h]
            else:
                piv[h] = cur
                rank += 1
                break
    return rank


def betti_gf2(maximal, maxdim=3):
    """Betti numbers b_0..b_maxdim over GF(2)."""
    by = _all_faces(maximal, maxdim + 1)
    idx = {d: {s: i for i, s in enumerate(v)} for d, v in by.items()}
    rank = {}
    for d in sorted(by):
        if d - 1 in by:
            cols = []
            for s in by[d]:
                bits = 0
                for i in range(len(s)):
                    bits ^= 1 << idx[d - 1][s[:i] + s[i + 1:]]
                cols.append(bits)
            rank[d] = _gf2_rank(cols)
        else:
            rank[d] = 0
    return tuple(len(by.get(d, [])) - rank.get(d, 0) - rank.get(d + 1, 0)
                 for d in range(maxdim + 1))


def euler_char(maximal, maxdim=4):
    by = _all_faces(maximal, maxdim)
    return sum((-1) ** d * len(v) for d, v in by.items())


# ======================================================================
#  stability scan and Morse pairing
# ======================================================================
def stability_scan(R, A, caps=(2, 3, 4, 5, 6, 8, 10, 12), nthick=4,
                   maxdim=3, hard_cap=13):
    """Betti vector versus shadow-size resolution.  The plateau is the signal.

    caps : the resolution parameter -- maximum cardinality of a cover set.
           Small caps under-resolve (spurious components); large caps
           over-connect (holes filled in).  A manifoldlike causet shows a
           plateau in between.
    """
    layer = layer_index(R, A, nmax=nthick + 1)
    out = {}
    for c in caps:
        mx = mrs_nerve_fast(R, A, layer, nthick, shadow_cap=c)
        big = max((len(s) for s in mx), default=0)
        if big > hard_cap:
            out[c] = ("too_big", big, len(mx))
            continue
        out[c] = (betti_gf2(mx, maxdim), big, len(mx))
    return out


def morse_count(beta_sequence):
    """Signed count of topology-change events from a foliation.

    beta_sequence : list of Betti tuples, one per antichain, in time order.
    A change in b_k by +1 is an index-k handle attachment; by -1 its dual.
    Returns sum_i (-1)^lambda_i  and the list of events.
    """
    events = []
    for (b0, b1) in zip(beta_sequence[:-1], beta_sequence[1:]):
        if b0 is None or b1 is None:
            continue
        for k, (u, v) in enumerate(zip(b0, b1)):
            if v > u:
                events += [(k, +1)] * (v - u)
            elif v < u:
                events += [(k, -1)] * (u - v)
    chi = sum((-1) ** k for k, s in events)
    return chi, events
