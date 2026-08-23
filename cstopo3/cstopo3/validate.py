"""
cstopo3.validate -- exact triangulations with known Betti numbers.

These test the homology engine, not the physics.  The essential ones are the
two manifolds the 3+1D demonstration must distinguish, S^3 and S^1 x S^2, and
the wedge S^1 v S^2 v S^3 which has the same Betti vector as S^1 x S^2 and is
therefore the negative control for the manifold condition.
"""
from __future__ import annotations
from itertools import combinations

__all__ = ["sphere_boundary_simplex", "octahedron_s2", "prism_over_complex",
           "s1_times_s2", "torus3", "wedge_s1_s2_s3", "KNOWN"]


def sphere_boundary_simplex(n):
    """S^{n-1} as the boundary of the n-simplex.  beta = (1,0,...,0,1)."""
    return [tuple(c) for c in combinations(range(n + 1), n)]


def octahedron_s2():
    """S^2 as the octahedron: 6 vertices, 8 triangles.  beta = (1,0,1)."""
    p = [0, 1, 2, 3, 4, 5]          # +x,-x,+y,-y,+z,-z
    tris = []
    for a in (0, 1):
        for b in (2, 3):
            for c in (4, 5):
                tris.append(tuple(sorted((p[a], p[b], p[c]))))
    return tris


def prism_over_complex(facets, nlevels, closed=True):
    """Triangulate K x S^1 (closed=True) or K x [0,1] (closed=False).

    Uses the standard staircase triangulation of a prism over a d-simplex:
    for sigma = (v_0 < ... < v_d) with copies at levels L and L+1, the prism
    splits into the d+1 simplices
        (v_0^L, ..., v_i^L, v_i^{L+1}, ..., v_d^{L+1}),   i = 0..d,
    which is a valid triangulation agreeing on shared faces because the vertex
    order is global.
    """
    out = []
    nl = nlevels if closed else nlevels - 1
    for L in range(nl):
        M = (L + 1) % nlevels if closed else L + 1
        for s in facets:
            s = tuple(sorted(s))
            d = len(s) - 1
            for i in range(d + 1):
                cell = tuple([(v, L) for v in s[:i + 1]] +
                             [(v, M) for v in s[i:]])
                # drop the degenerate repeat when L == M (nlevels == 1)
                if len(set(cell)) == len(s) + 1:
                    out.append(tuple(sorted(set(cell))))
    return out


def s1_times_s2(nlevels=3):
    """S^1 x S^2 as a prism over the octahedron, ends identified.

    beta = (1,1,1,1).  This is the spatial slice created by the index-one
    surgery of the wormhole sector.
    """
    return prism_over_complex(octahedron_s2(), nlevels, closed=True)


def torus3(n=3):
    """T^3 = (S^1)^3 by iterated prism.  beta = (1,3,3,1)."""
    circle = [(i, (i + 1) % n) for i in range(n)]        # S^1, n>=3
    t2 = prism_over_complex(circle, n, closed=True)      # T^2
    return prism_over_complex(t2, n, closed=True)        # T^3


def wedge_s1_s2_s3():
    """S^1 v S^2 v S^3, beta = (1,1,1,1): same Betti vector as S^1 x S^2.

    Built as boundary-of-4-simplex (S^3 on vertices 0..4) with one extra vertex
    5 joined to an edge (creating a 1-cycle) and one extra vertex 6 coned over
    a triangle boundary (creating a 2-cycle).
    """
    f = [tuple(c) for c in combinations(range(5), 4)]     # S^3
    f += [(0, 5), (1, 5)]                                # S^1 wedge
    f += [(0, 1, 6), (1, 2, 6), (0, 2, 6)]               # S^2 wedge
    return f


KNOWN = {
    "S^2 (octahedron)":      (octahedron_s2(),               (1, 0, 1)),
    "S^3 (bdy 4-simplex)":   (sphere_boundary_simplex(4),    (1, 0, 0, 1)),
    "S^4 (bdy 5-simplex)":   (sphere_boundary_simplex(5),    (1, 0, 0, 0, 1)),
    "S^1 x S^2":             (s1_times_s2(3),                (1, 1, 1, 1)),
    "S^1 x S^2 (4 levels)":  (s1_times_s2(4),                (1, 1, 1, 1)),
    "T^3":                   (torus3(3),                     (1, 3, 3, 1)),
    "S^1 v S^2 v S^3":       (wedge_s1_s2_s3(),              (1, 1, 1, 1)),
}
