"""
Delta-beta search
=================

Search finite causal sets for the wormhole spatial transition

    S^3  ->  S^1 x S^2  ->  S^3
    (1,0,0,1) -> (1,1,1,1) -> (1,0,0,1)

detected through nerves of thickened antichains (Major-Rideout-Surya).

IMPORTANT CAVEAT (see Sec. III.C of the paper):
    beta = (1,1,1,1) does NOT identify S^1 x S^2.  The wedge S^1 v S^2 v S^3
    has the same Betti vector and is distinguished only by the cup product
    a \smile b.  A hit from this search is therefore a CANDIDATE, not a
    wormhole; the manifold / cup-product condition must be checked separately
    (see `is_pseudomanifold` below for a cheap necessary test).

Usage:
    python delta_beta_search.py --nmax 18 --trials 200000
"""

import numpy as np, itertools, argparse
from collections import Counter
from homology import betti_and_chi


# ----------------------------------------------------------------------
# causal sets
# ----------------------------------------------------------------------
def transitive_closure(R):
    n = R.shape[0]
    for k in range(n):
        R |= np.outer(R[:, k], R[k, :])
    return np.triu(R, 1)


def random_causet(n, p, rng):
    return transitive_closure(np.triu(rng.random((n, n)) < p, 1))


def layered_causet(n_layers, width, p, rng):
    """A causet with an explicit layer structure, closer to a sprinkling
       of a thin slab than transitive percolation is."""
    n = n_layers * width
    R = np.zeros((n, n), bool)
    for a in range(n_layers - 1):
        for b in range(a + 1, n_layers):
            for i in range(width):
                for j in range(width):
                    if rng.random() < p / (b - a):
                        R[a * width + i, b * width + j] = True
    return transitive_closure(R)


# ----------------------------------------------------------------------
# MRS-style nerve of a thickened antichain
# ----------------------------------------------------------------------
def maximal_antichains(R, min_size=4):
    n = R.shape[0]
    comp = R | R.T
    adj = [set(np.nonzero(~comp[i])[0]) - {i} for i in range(n)]
    out = []

    def bk(Rset, P, X):
        if not P and not X:
            if len(Rset) >= min_size:
                out.append(sorted(Rset))
            return
        if not P:
            return
        piv = max(P | X, key=lambda v: len(adj[v] & P))
        for v in list(P - adj[piv]):
            bk(Rset | {v}, P & adj[v], X & adj[v])
            P = P - {v}
            X = X | {v}

    bk(set(), set(range(n)), set())
    return out


def nerve_maximal_simplices(R, A, thickening):
    """Cover of the antichain A by past-shadows of elements lying within
       `thickening` links to the future.  Returns the maximal simplices."""
    n = R.shape[0]
    Aarr = np.array(A)
    # layer index above A
    layer = np.full(n, -1, int)
    layer[Aarr] = 0
    order = np.argsort([R[:, x].sum() for x in range(n)])
    for x in order:
        if layer[x] == 0:
            continue
        preds = np.nonzero(R[:, x])[0]
        best = -1
        for q in preds:
            if layer[q] >= 0 and layer[q] + 1 > best:
                best = layer[q] + 1
        layer[x] = best
    slab = np.nonzero((layer >= 1) & (layer <= thickening))[0]

    shadows = set()
    for y in slab:
        s = tuple(np.nonzero(R[Aarr, y])[0])
        if s:
            shadows.add(s)
    # keep only maximal shadows
    mx = [s for s in shadows if not any(set(s) < set(t) for t in shadows)]
    covered = set()
    for s in mx:
        covered.update(s)
    for i in range(len(A)):
        if i not in covered:
            mx.append((i,))
    return mx


def betti_vector(maximal, kmax=3):
    b, chi, f = betti_and_chi(maximal)
    return tuple(b.get(k, 0) for k in range(kmax + 1))


def is_pseudomanifold(maximal, d=3):
    """cheap necessary condition for a d-manifold: every (d-1)-face of a
       d-simplex lies in exactly two d-simplices."""
    top = [s for s in maximal if len(s) == d + 1]
    if not top:
        return False
    cnt = Counter()
    for s in top:
        for f in itertools.combinations(sorted(s), d):
            cnt[f] += 1
    return all(v == 2 for v in cnt.values())


# ----------------------------------------------------------------------
# search
# ----------------------------------------------------------------------
TARGET_S3 = (1, 0, 0, 1)
TARGET_WH = (1, 1, 1, 1)


def scan_causet(R, thickenings=(1, 2, 3)):
    """Return the set of (antichain, thickening, betti) triples observed."""
    hits = []
    for A in maximal_antichains(R):
        for th in thickenings:
            mx = nerve_maximal_simplices(R, A, th)
            if max((len(s) for s in mx), default=0) > 9:
                continue                      # complex too large to homologise
            try:
                bv = betti_vector(mx)
            except Exception:
                continue
            hits.append((tuple(A), th, bv, mx))
    return hits


def search(nmax, trials, seed=0, verbose=True):
    rng = np.random.default_rng(seed)
    best = None
    stats = Counter()
    for t in range(trials):
        n = rng.integers(8, nmax + 1)
        p = rng.uniform(0.15, 0.5)
        R = random_causet(int(n), float(p), rng)
        for A, th, bv, mx in scan_causet(R):
            stats[bv] += 1
            if bv == TARGET_WH:
                pm = is_pseudomanifold(mx)
                if best is None or n < best[0]:
                    best = (int(n), tuple(A), th, bv, pm)
                    if verbose:
                        print(f"  candidate: n={n}, |A|={len(A)}, thickening={th}, "
                              f"beta={bv}, pseudomanifold={pm}")
    return best, stats


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=18)
    ap.add_argument("--trials", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    print("=" * 74)
    print(" Delta-beta search for the wormhole spatial transition")
    print("=" * 74)
    print(f" target intermediate beta = {TARGET_WH}   (S^3 endpoints: {TARGET_S3})")
    print(" NOTE: beta=(1,1,1,1) is necessary, NOT sufficient -- wedges qualify.\n")

    best, stats = search(args.nmax, args.trials, args.seed)

    print("\n most frequent Betti vectors observed:")
    for bv, c in stats.most_common(8):
        print(f"    {bv}: {c}")
    print(f"\n best candidate: {best}")
    if best is None:
        print(" => no realisation found in this range; consistent with the")
        print("    estimate n_WH >~ 17 (Betti signature) or >~ 40 (genuine S^1xS^2).")
