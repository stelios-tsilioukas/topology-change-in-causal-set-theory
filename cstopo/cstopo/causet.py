"""
cstopo.causet -- sprinkling and causal relations.

Dependency-light: numpy only.

The causal matrix is stored as a packed bit array (np.uint64) to keep memory at
N^2/8 bytes.  Relations are built with a light-cone time window so the cost is
O(N * <neighbours>) rather than O(N^2) for elongated regions.
"""
from __future__ import annotations
import numpy as np

__all__ = ["Sprinkling", "CausalMatrix", "sprinkle_box", "build_causal_flat"]


# ----------------------------------------------------------------------
# packed boolean matrix
# ----------------------------------------------------------------------
class CausalMatrix:
    """Strict causal relation R[i,j] = True iff element i precedes element j.

    Stored row-wise as packed bits.  Only the upper triangle in time order is
    populated (elements are sorted by time on construction), so R is a strict
    partial order by construction and needs no transitive closure when built
    from an embedding.
    """

    __slots__ = ("n", "_bits", "_words")

    def __init__(self, n: int):
        self.n = n
        self._words = (n + 63) // 64
        self._bits = np.zeros((n, self._words), dtype=np.uint64)

    # --- element access -------------------------------------------------
    def set_row(self, i: int, js: np.ndarray) -> None:
        if js.size == 0:
            return
        w = (js >> 6).astype(np.int64)
        b = (js & 63).astype(np.uint64)
        np.bitwise_or.at(self._bits[i], w, np.uint64(1) << b)

    def row(self, i: int) -> np.ndarray:
        """Indices of elements to the strict future of i."""
        return np.flatnonzero(np.unpackbits(
            self._bits[i].view(np.uint8)).reshape(-1, 8)[:, ::-1].ravel()[: self.n])

    def row_mask(self, i: int) -> np.ndarray:
        return np.unpackbits(
            self._bits[i].view(np.uint8)).reshape(-1, 8)[:, ::-1].ravel()[: self.n].astype(bool)

    def rows_and(self, i: int, j: int) -> np.ndarray:
        """Elements in the future of both i and j."""
        return np.flatnonzero(np.unpackbits(
            (self._bits[i] & self._bits[j]).view(np.uint8)
        ).reshape(-1, 8)[:, ::-1].ravel()[: self.n])

    def popcount_row(self, i: int) -> int:
        v = self._bits[i]
        return int(sum(bin(int(x)).count("1") for x in v))

    @property
    def nbytes(self) -> int:
        return self._bits.nbytes

    def to_dense(self) -> np.ndarray:
        """Dense boolean matrix.  Only for small n -- diagnostics/tests."""
        out = np.zeros((self.n, self.n), dtype=bool)
        for i in range(self.n):
            out[i] = self.row_mask(i)
        return out


# ----------------------------------------------------------------------
# sprinkling
# ----------------------------------------------------------------------
class Sprinkling:
    """A Poisson sprinkling: coordinates plus the causal matrix.

    Attributes
    ----------
    t : (N,) proper/coordinate time, sorted ascending
    x : (N, d-1) spatial coordinates
    R : CausalMatrix
    """

    __slots__ = ("t", "x", "R", "meta")

    def __init__(self, t, x, R, meta=None):
        self.t, self.x, self.R = t, x, R
        self.meta = meta or {}

    @property
    def N(self):
        return len(self.t)

    @property
    def dim(self):
        return self.x.shape[1] + 1


def sprinkle_box(N, T, L, d=4, seed=0):
    """Uniform sprinkling into a box: t in [0,T), x in [-L/2, L/2)^(d-1)."""
    rng = np.random.default_rng(seed)
    t = np.sort(rng.random(N) * T)
    x = rng.random((N, d - 1)) * L - L / 2
    return t, x


def build_causal_flat(t, x, chunk=4096):
    """Causal matrix for flat space: i < j iff (t_j - t_i) > |x_j - x_i|.

    Elements must be time-sorted.  Uses a time window: for element i we need
    only consider j with t_j - t_i <= max possible spatial separation, but in a
    box of side L that is L*sqrt(d-1), so we window on that.
    """
    N = len(t)
    R = CausalMatrix(N)
    # maximum spatial separation in the box bounds the useful time window
    span = np.sqrt(((x.max(0) - x.min(0)) ** 2).sum())
    for i0 in range(0, N, chunk):
        i1 = min(i0 + chunk, N)
        # candidate futures: j > i with t_j - t_i <= span is NOT a restriction
        # (any j with larger dt is automatically related); handle both parts.
        for i in range(i0, i1):
            jhi = np.searchsorted(t, t[i] + span, side="right")
            j0 = i + 1
            if j0 >= N:
                continue
            # part 1: j in (i, jhi) -- test explicitly
            if jhi > j0:
                dt = t[j0:jhi] - t[i]
                dx2 = ((x[j0:jhi] - x[i]) ** 2).sum(1)
                rel = np.flatnonzero(dt * dt > dx2) + j0
            else:
                rel = np.empty(0, dtype=np.int64)
            # part 2: j >= jhi -- always related (dt > span >= |dx|)
            if jhi < N:
                rel = np.concatenate([rel, np.arange(jhi, N)])
            R.set_row(i, rel)
    return R


def ordering_fraction(R):
    """Fraction of pairs that are related -- a basic geometric diagnostic."""
    n = R.n
    tot = sum(R.popcount_row(i) for i in range(n))
    return 2.0 * tot / (n * (n - 1)) / 2.0 * 2.0 if n > 1 else 0.0
