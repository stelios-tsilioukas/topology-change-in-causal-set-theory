"""
cstopo3.geom3 -- 3+1D ultrastatic spacetimes with three-dimensional slices.

The causal relation is  i < j  iff  t_j - t_i > d_spatial(x_i, x_j), exact for
ds^2 = -dt^2 + h_ij dx^i dx^j.  Coordinates are used ONCE, to generate the
relation; every topological step downstream sees only `precedes`.

Two static slices with *exact* closed-form geodesic distance:

    S3            round 3-sphere,      beta = (1,0,0,1)
    S1xS2         Riemannian product,  beta = (1,1,1,1)

and one cobordism between them, built as a warped product over the polar angle,
which is the index-one surgery of the wormhole sector:

    ds^2_slice = dpsi^2 + r_lambda(psi)^2 dOmega^2,   r_lambda = sqrt(sin^2 psi + lambda^2)

    lambda = 0, ends are poles                     ->  S^3
    lambda > 0, ends are 2-spheres, glued to each other ->  S^1 x S^2

The surgery is therefore a single parameter: opening the two poles into
2-spheres of radius lambda and identifying them.  Distances on a warped product
have no closed form; we use the product approximation with the mean warp
factor, as in the 2+1D construction.  `S3Warped` exists so that this
approximation can be validated against the exact `S3` through the identical
pipeline.
"""
from __future__ import annotations
import numpy as np

# np.trapezoid is the NumPy >= 2.0 spelling; np.trapz is the < 2.0 one, removed
# in 2.0.  Bind whichever exists so the package runs on both generations.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz
__all__ = ["S3", "S1xS2", "S3Warped", "Surgery", "Spacetime"]


def _unit_sphere(n, dim, rng):
    """n uniform points on the unit sphere in R^dim."""
    v = rng.standard_normal((n, dim))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def _ang(u, v):
    """Angle between unit vectors; u is (dim,) or (n,dim), v is (m,dim)."""
    return np.arccos(np.clip(u @ v.T, -1.0, 1.0))


# ======================================================================
#  static slices, exact distances
# ======================================================================
class S3:
    """Round 3-sphere of radius R.  Exact geodesic distance."""
    name = "S3"
    betti = (1, 0, 0, 1)

    def __init__(self, R=1.0):
        self.R = R

    @property
    def volume(self):
        return 2 * np.pi ** 2 * self.R ** 3

    def sample(self, n, rng):
        return _unit_sphere(n, 4, rng)

    def dist(self, x, i, js):
        return self.R * _ang(x[i], x[js])


class S1xS2:
    """Riemannian product S^1(R1) x S^2(R2).  Exact: d^2 = d_1^2 + d_2^2."""
    name = "S1xS2"
    betti = (1, 1, 1, 1)

    def __init__(self, R1=1.0, R2=1.0):
        self.R1, self.R2 = R1, R2

    @property
    def volume(self):
        return (2 * np.pi * self.R1) * (4 * np.pi * self.R2 ** 2)

    def sample(self, n, rng):
        phi = rng.random(n) * 2 * np.pi
        om = _unit_sphere(n, 3, rng)
        return np.column_stack([phi, om])          # (n,4): phi, omega

    def dist(self, x, i, js):
        dphi = np.abs(x[js, 0] - x[i, 0])
        dphi = np.minimum(dphi, 2 * np.pi - dphi)
        dth = _ang(x[i, 1:], x[js, 1:])
        return np.sqrt((self.R1 * dphi) ** 2 + (self.R2 * dth) ** 2)


# ======================================================================
#  warped product: the surgery family
# ======================================================================
class _Warped:
    """ds^2 = dpsi^2 + r(psi)^2 dOmega^2 for psi in [psi_lo, psi_hi].

    `glue` identifies the two ends, making psi a circle of circumference
    psi_hi - psi_lo.  Distance uses the product approximation with the mean warp
    factor: exact where r is constant, and correct in the pinched limit.
    """
    name = "warped"

    def __init__(self, R=1.0, lam=0.0, neck=0.0, glue=False):
        self.R, self.lam, self.neck, self.glue = R, lam, neck, glue
        if glue:
            self.psi_lo, self.psi_hi = 0.0, 2 * neck + np.pi
        else:
            self.psi_lo, self.psi_hi = neck, neck + np.pi

    def r(self, psi):
        """lam-thickened S^3 profile with straight tube sections of length
        `neck` at each end.  For psi below neck or above neck+pi the sine
        argument saturates and r = lam, a genuine cylindrical tube."""
        s = np.sin(np.clip(np.asarray(psi, float) - self.neck, 0.0, np.pi))
        return np.sqrt(s ** 2 + self.lam ** 2)

    @property
    def circumference(self):
        return self.R * (self.psi_hi - self.psi_lo) if self.glue else np.inf

    @property
    def volume(self):
        psi = np.linspace(self.psi_lo, self.psi_hi, 40001)
        return 4 * np.pi * self.R ** 3 * _trapezoid(self.r(psi) ** 2, psi)

    def sample(self, n, rng):
        grid = np.linspace(self.psi_lo, self.psi_hi, 40001)
        cdf = np.cumsum(self.r(grid) ** 2); cdf /= cdf[-1]
        psi = np.interp(rng.random(n), cdf, grid)
        return np.column_stack([psi, _unit_sphere(n, 3, rng)])

    def dist(self, x, i, js):
        dpsi = np.abs(x[js, 0] - x[i, 0])
        if self.glue:
            span = self.psi_hi - self.psi_lo
            dpsi = np.minimum(dpsi, span - dpsi)
        rbar = 0.5 * (self.r(x[i, 0]) + self.r(x[js, 0]))
        dth = _ang(x[i, 1:], x[js, 1:])
        return self.R * np.sqrt(dpsi ** 2 + (rbar * dth) ** 2)


class S3Warped(_Warped):
    """S^3 in warped coordinates, occupying psi in [neck, neck+pi] so that it
    sits inside the coordinate range of the post-surgery slice.  Control for the
    product approximation: must reproduce the exact S3 result."""
    name = "S3warped"
    betti = (1, 0, 0, 1)

    def __init__(self, R=1.0, neck=0.0):
        super().__init__(R=R, lam=0.0, neck=neck, glue=False)


class Surgery(_Warped):
    """Post-surgery slice: poles opened to throat radius lam, a tube of length
    `neck` inserted at each end, and the ends identified.  Topologically
    S^1 x S^2; the 1-cycle has length 2*neck + pi, which must exceed the cover
    radius or it is filled in."""
    name = "surgery"
    betti = (1, 1, 1, 1)

    def __init__(self, R=1.0, lam=0.35, neck=0.0):
        super().__init__(R=R, lam=lam, neck=neck, glue=True)


# ======================================================================
#  spacetime: sprinkling and the causal relation
# ======================================================================
class Spacetime:
    """Ultrastatic slab Sigma x [0,T], or a cobordism at t_split.

    slice_before / slice_after must share a coordinate convention.  For a
    static run pass only `slice_before`.
    """

    def __init__(self, slice_before, N, T=1.0, seed=0,
                 slice_after=None, t_split=None):
        self.sb = slice_before
        self.sa = slice_after if slice_after is not None else slice_before
        self.t_split = t_split if t_split is not None else np.inf
        self.T = T
        rng = np.random.default_rng(seed)
        # A Poisson sprinkling has constant density per unit FOUR-volume, so
        # when the slice volume jumps at t_split the number of elements per unit
        # time must jump with it.
        if np.isfinite(self.t_split) and 0 < self.t_split < T:
            Vb = self.sb.volume * self.t_split
            Va = self.sa.volume * (T - self.t_split)
            Nb = int(round(N * Vb / (Vb + Va)))
            tb = np.sort(rng.random(Nb) * self.t_split)
            ta = np.sort(self.t_split + rng.random(N - Nb) * (T - self.t_split))
            t = np.concatenate([tb, ta])
        else:
            t = np.sort(rng.random(N) * T)
        self.t = t
        self.N = N
        pre = t < self.t_split
        x = np.empty((N, 4))
        if pre.any():
            x[pre] = self.sb.sample(int(pre.sum()), rng)
        if (~pre).any():
            x[~pre] = self.sa.sample(int((~pre).sum()), rng)
        self.x = x
        self.pre = pre
        # discreteness scale from the four-volume actually used
        vol = self.sb.volume * min(T, self.t_split if np.isfinite(self.t_split) else T)
        if np.isfinite(self.t_split):
            vol += self.sa.volume * max(0.0, T - self.t_split)
        self.V4 = vol
        self.ell = (vol / N) ** 0.25

    def sdist(self, i, js):
        """Spatial distance, evaluated on the slice of the later element.

        A causal path runs forward in time, so it is the later geometry that
        controls whether a shortcut exists.  For a static run this is exact.
        """
        js = np.asarray(js)
        later = self.sa if self.t[i] >= self.t_split else self.sb
        if np.isfinite(self.t_split):
            out = np.empty(js.size)
            m = self.t[js] >= self.t_split
            if m.any():
                out[m] = self.sa.dist(self.x, i, js[m])
            if (~m).any():
                out[~m] = self.sb.dist(self.x, i, js[~m])
            return out
        return later.dist(self.x, i, js)

    def precedes(self, i, js):
        """Boolean: does element i causally precede each of js?"""
        js = np.asarray(js)
        return (self.t[js] - self.t[i]) > self.sdist(i, js)
