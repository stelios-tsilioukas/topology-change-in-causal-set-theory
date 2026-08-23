"""
cstopo.geometry -- spacetimes to sprinkle into.

Each generator returns (t, x, R, meta) with R a CausalMatrix built from the
*causal structure of that geometry*, after which the embedding may be discarded:
everything downstream uses only the order.

Two geometries are provided:

  trousers_2p1   a 2+1 dimensional topology-changing cobordism whose spatial
                 slices go  one circle -> two circles.  The analogue of
                 S^3 -> S^1 x S^2, but with b_1 directly meaningful.

  gs_wormhole    a 3+1 dimensional Giddings-Strominger wormhole,
                 ds^2 = dtau^2 + a(tau)^2 dOmega_3^2,  a^4 = a0^4 + tau^2...
                 (see below), whose spatial slices go S^3 -> S^1 x S^2 -> S^3.
"""
from __future__ import annotations
import numpy as np
from .causet import CausalMatrix

__all__ = ["trousers_2p1", "gs_wormhole", "cylinder_2p1"]


# ======================================================================
#  generic builder
# ======================================================================
def _build_from_metric(t, x, dist_fn, chunk_report=None):
    """Causal matrix from a spatial distance function on time-sorted points.

    i precedes j iff  t_j - t_i > d_spatial(i, j).  This is exact for a
    static/ultrastatic metric ds^2 = -dt^2 + h_ij dx^i dx^j and is the standard
    approximation for slowly varying a(t) provided the sprinkling is dense
    compared with the scale on which a varies.
    """
    N = len(t)
    R = CausalMatrix(N)
    for i in range(N):
        j0 = i + 1
        if j0 >= N:
            continue
        d = dist_fn(i, np.arange(j0, N))
        dt = t[j0:] - t[i]
        rel = np.flatnonzero(dt > d) + j0
        R.set_row(i, rel)
    return R


# ======================================================================
#  2+1D trousers
# ======================================================================
def _trousers_radius(tt, t_split, r_waist, r_leg, sep):
    """Half-separation of the two leg centres and the tube radius at time tt.

    For tt < t_split the slice is a single circle of radius r_waist.
    For tt > t_split it is two circles of radius r_leg whose centres are at
    +/- s(tt), with s growing from 0.
    """
    if tt <= t_split:
        return 0.0, r_waist
    frac = (tt - t_split)
    return sep * frac, r_leg


def trousers_2p1(N, T=2.0, t_split=1.0, r_waist=1.0, r_leg=0.6,
                 sep=1.2, seed=0):
    """Sprinkle a 2+1D trousers.

    The spatial slice at time t is a 1-manifold (circle, or two circles)
    embedded in the plane; the spacetime metric is taken ultrastatic,
    ds^2 = -dt^2 + dl^2 with dl the arclength along that 1-manifold.
    Causal relation: t_j - t_i > geodesic distance along the slice.

    Returns (t, x, R, meta).  x holds the embedding (for diagnostics only).
    """
    rng = np.random.default_rng(seed)
    t = np.sort(rng.random(N) * T)

    # place each element on the slice appropriate to its time
    theta = rng.random(N) * 2 * np.pi
    comp = np.zeros(N, dtype=np.int8)          # which leg (0 = waist / left)
    cx = np.zeros(N)
    rad = np.zeros(N)
    for k in range(N):
        s, r = _trousers_radius(t[k], t_split, r_waist, r_leg, sep)
        rad[k] = r
        if s == 0.0:
            cx[k] = 0.0
        else:
            comp[k] = rng.integers(0, 2)
            cx[k] = -s if comp[k] == 0 else +s
    x = np.column_stack([cx + rad * np.cos(theta), rad * np.sin(theta)])

    def dist(i, js):
        """Geodesic distance along the slice.

        Same component: arclength r*|dtheta| (shorter way round).
        Different components after the split: infinite (spatially disconnected).
        """
        d = np.empty(len(js))
        dth = np.abs(theta[js] - theta[i])
        dth = np.minimum(dth, 2 * np.pi - dth)
        rmean = 0.5 * (rad[i] + rad[js])
        d[:] = rmean * dth
        # disconnected if both are past the split and on different legs
        split_i = t[i] > t_split
        split_j = t[js] > t_split
        if split_i:
            bad = split_j & (comp[js] != comp[i])
            d[bad] = np.inf
        return d

    R = _build_from_metric(t, x, dist)
    meta = dict(kind="trousers_2p1", T=T, t_split=t_split, r_waist=r_waist,
                r_leg=r_leg, sep=sep, N=N, seed=seed)
    return t, x, R, meta


def cylinder_2p1(N, T=2.0, r=1.0, seed=0):
    """Control geometry: a 2+1D cylinder, spatial slice always one circle."""
    rng = np.random.default_rng(seed)
    t = np.sort(rng.random(N) * T)
    theta = rng.random(N) * 2 * np.pi
    x = np.column_stack([r * np.cos(theta), r * np.sin(theta)])

    def dist(i, js):
        dth = np.abs(theta[js] - theta[i])
        dth = np.minimum(dth, 2 * np.pi - dth)
        return r * dth

    R = _build_from_metric(t, x, dist)
    return t, x, R, dict(kind="cylinder_2p1", T=T, r=r, N=N, seed=seed)


# ======================================================================
#  3+1D Giddings-Strominger wormhole
# ======================================================================
def gs_wormhole(N, a0=1.0, tau_max=3.0, seed=0):
    """Sprinkle the Euclidean GS wormhole, ds^2 = dtau^2 + a(tau)^2 dOmega_3^2.

    From (a')^2 = 1 - a0^4/a^4 the throat is at a = a0.  We parametrise by
    tau in [-tau_max, tau_max] and solve a(tau) by quadrature.

    Spatial slices are 3-spheres of radius a(tau) except near the throat, where
    the geometry is a neck: the slice topology in the *Lorentzian* continuation
    is S^3 away from the throat and the handle is the neck itself.  For the
    purposes of the homology test we sprinkle the neck region and read the
    slice topology directly.
    """
    rng = np.random.default_rng(seed)

    # a(tau) by quadrature: dtau = da / sqrt(1 - a0^4/a^4)
    agrid = np.linspace(a0, a0 * 8, 4000)
    integ = 1.0 / np.sqrt(np.clip(1 - (a0 / agrid) ** 4, 1e-12, None))
    taugrid = np.concatenate([[0.0], np.cumsum(0.5 * (integ[1:] + integ[:-1]) *
                                               np.diff(agrid))])
    tau_max = min(tau_max, taugrid[-1])

    # sample tau with density proportional to the proper 4-volume a^3
    tt = np.linspace(-tau_max, tau_max, 4000)
    aa = np.interp(np.abs(tt), taugrid, agrid)
    w = aa ** 3
    cdf = np.cumsum(w); cdf /= cdf[-1]
    tau = np.interp(rng.random(N), cdf, tt)
    tau = np.sort(tau)
    a = np.interp(np.abs(tau), taugrid, agrid)

    # uniform points on S^3
    v = rng.normal(size=(N, 4))
    v /= np.linalg.norm(v, axis=1, keepdims=True)
    x = v * a[:, None]

    def dist(i, js):
        # geodesic distance on the S^3 of the mean radius
        dot = np.clip((v[js] * v[i]).sum(1), -1.0, 1.0)
        amean = 0.5 * (a[i] + a[js])
        return amean * np.arccos(dot)

    R = _build_from_metric(tau, x, dist)
    meta = dict(kind="gs_wormhole", a0=a0, tau_max=float(tau_max), N=N,
                seed=seed, a_of_tau=(taugrid, agrid))
    return tau, x, R, meta
