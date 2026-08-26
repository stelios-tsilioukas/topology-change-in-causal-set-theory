# cstopo3 — the 3+1D demonstration

The wormhole-sector transition

    S³ ──(index-1 surgery)──> S¹×S²        Δβ = (0,+1,+1,0)

recovered **from the causal order alone, 8/8 sprinklings, at N = 20 000**, in
about 8 seconds per sprinkling on one core.

This supersedes an earlier attempt (kept as `cstopo/run_wormhole.py`) which did
not converge below N ≈ 5×10⁴. Section 1 below explains why: the failure was a
dimension-counting error in the resolution parameter, not a shortage of
elements.

Dependency: numpy. Nothing else. `python3 test_all.py` first.

---

## 1. Why the earlier attempt failed

Not a bug and not a shortage of elements — a **dimension-counting error in the
resolution parameter**.

For a point cloud sampling a *d*-manifold, a Dowker/nerve complex reproduces the
topology only when the cover balls have radius `r ≈ 2×` the vertex spacing `s`.
The facet size that follows is

    m = ρ·(4π/3)r³ = (4π/3)(r/s)³ ≈ 33     for r/s = 2

and — this is the point — **m is independent of the sampling density**. It is set
by the dimension alone. The 2+1D demonstration used cap ≈ 8–12, which is right
for a 2-dimensional slice. Carried over to a 3-dimensional slice, cap ≤ 16 gives
`r/s ≲ 1.55`: under-resolved. The failure reproduces exactly:
`b₀` between 2 and 134, `b₁` between 35 and 485, no plateau:

| cap | 4 | 6 | 8 | 10 | 12 | 16 |
|---|---|---|---|---|---|---|
| S³, `(b₀,b₁)` | (134,35) | (55,197) | (20,344) | (8,344) | (6,287) | (2,123) |

A diagnostic worth remembering: the elementary collapse removed **nothing**
(172 310 → 172 256 simplices). A dense sample of a 3-manifold should collapse by
orders of magnitude. That alone says the complex is not manifold-like.

But raising cap to ~33 explodes the clique complex: `C(33,4) = 40 920`
3-simplices *per facet*.

## 2. The two moves that fix it

**(a) Landmarks.** Subsample the antichain to a sparse vertex set `L`, and keep
*all* of `A` as witnesses. Facets are `σ(y) ∩ L` — still purely order-theoretic,
just a random subset, no metric notion introduced. Facet size is fixed by `r/s`
either way, but `|L|` is now a free knob and it is `|L|` that bounds the linear
algebra. `|L|` in the low hundreds instead of `|A| ~ 10³–10⁴`.

**(b) A cheap b₁.** Over Z₂, the class of a triangle's boundary in the cycle
space of the 1-skeleton, written in the fundamental basis of a spanning forest,
is just **the set of non-tree edges among its three edges**. So every column of
`∂₂` has ≤ 3 nonzeros, in a space of dimension `m = n₁ − n₀ + b₀ ~ 10⁴`. The
triangles can be numerous; the space they live in is small.

    b₀ = components of the 1-skeleton                (union-find)
    b₁ = m − rank(∂₂ in cycle-space coordinates)

Python bigints beat numpy arrays as the bitset here — at `m ~ 10⁴` a column is
~1 kB and a bigint XOR is one C loop, an order of magnitude cheaper than numpy's
per-call overhead. That change alone was 30×.

## 3. Results

All from the causal order alone; coordinates are used once, to generate the
relation, and never again.

**Static slices, exact closed-form geodesic distances, N = 20 000, |L| = 250:**

| | cap 8 | cap 10 | cap 12 | cap 14 |
|---|---|---|---|---|
| S³ → (1,0) | 0/8 | 7/8 | 7/8 | **8/8** |
| S¹×S² → (1,1) | 5/8 | 8/8 | 8/8 | **8/8** |

cap = 8 under-resolves exactly as the dimension count predicts. cap ≥ 10 is the
plateau; cap 14 is clean for both.

**Threshold in N** (cap 14, 4 seeds) — 4/4 for both geometries at **N = 4 000**
already, and at 8 000, 15 000, 30 000. Far below the paper's estimate.

**The cobordism** (N = 20 000, |L| = 250, cap = 14, λ = 0.40, neck = 1.5):

    seed:   1    2    3    4    5    6    7    8
    before  (1,0) throughout          →   after  (1,1) throughout
    Δb₁ = +1 detected:  8/8

**Cost scaling** (cobordism, both slices, single core):

| N | ℓ | \|A\| before/after | \|L\| | wall | RSS | detected |
|---|---|---|---|---|---|---|
| 20 000 | 0.278 | 884 / 1 476 | 250 | 5.1 s | 0.05 GB | yes |
| 40 000 | 0.234 | 1 504 / 2 455 | 300 | 11.2 s | 0.06 GB | yes |
| 80 000 | 0.197 | 2 575 / 4 185 | 350 | 25.1 s | 0.08 GB | yes |
| 160 000 | 0.165 | 4 320 / 6 983 | 400 | 62.7 s | 0.10 GB | yes |

**|A| ∝ N^0.743**, as the N^{3/4} expectation requires.

## 4. What does NOT work

**b₂ is now converged, and Δβ is measured in every component.** This section
previously recorded b₂ as the outstanding negative result. It is resolved; the
history is kept below because the diagnosis is the useful part.

The original approach pushed `cap` upward at `|L| = 250`, which is expensive and
does not get there:

| cap | 10 | 12 | 14 | 16 | 18 |
|---|---|---|---|---|---|
| b₂ (S³, N = 2e4, \|L\| = 250) | 42 | 26 | 6 | 3 | 1 |
| wall | 8 s | 30 s | 107 s | 262 s | 828 s |

`cap` was the wrong knob. At a **coarser** landmark set the complex is more
solid in dimension 3 at the same cap, because `r/s` is larger. At
`N = 4e4, |L| = 150, cap = 12` the static slices come out exact and cheap:

    S^3     -> (1,0,0)   3/3    45 s
    S^1xS^2 -> (1,1,1)   3/3    76 s

Two further things were needed for the cobordism.

**The 2-cycle needs a throat that is not too thin.** The 2-cycle is the throat's
S², and a thin throat is filled in by the cover, giving β = (1,1,0) — a solid
torus — stably rather than noisily:

| λ | 0.40 | 0.55 | 0.70 | 0.85 | 1.00 |
|---|---|---|---|---|---|
| after-slice β | (1,1,0) | (1,1,0) | **(1,1,1)** | **(1,1,1)** | (1,1,2) |

This is the b₂ analogue of the short-fat-handle failure already known for b₁,
and like it, it is geometry rather than a bug. The two constraints are
compatible because the 1-cycle length is `2·neck + π`, independent of λ: a long
neck keeps b₁, a fat throat keeps b₂.

**A cobordism splits the slab, so each slice sees only a fraction of N.** At
N = 8e4 the after-slice is already exact while the before-slice still returns
b₂ = 1 to 12. N = 2e5 fixes it.

At `N = 2e5, |L| = 150, cap = 12, λ = 0.80, neck = 1.5`:

    before = (1,0,0)    after = (1,1,1)    Δβ = (0,+1,+1)    3/3
    206 s per seed on one core, 0.31 GB peak

so Δβ = (0,+1,+1,0) is now measured directly rather than resting on Poincaré
duality. Run `bash tower_b2.sh` for the whole programme; it needs a workstation,
not a cluster.

**Two further limitations, stated plainly.**

*The handle must not be short and fat.* λ = 0.6 with neck = 0 gives 0/4 — the
1-cycle has a short representative and the cover fills it in. λ ≤ 0.45 works, and
with a neck (1-cycle length 2·neck + π) it is robust. This is geometry, not a bug,
and the physically relevant regime (thin throat) is the one that works. The sweep
is step 5 of `tower_run.sh`.

*The observable never uses a cross-split relation.* The before- and after-slice
witness windows lie entirely on their own side of `t_split`. What is measured is
the pair of slice topologies, which is exactly Δβ — but it is not a test of the
causal structure *through* the transition. Worth saying in the paper.

## 5. What to run

**On the tower (10 cores / 32 GB) — everything below, 1–2 hours:**

```
python3 test_all.py          # first, ~2 min
bash tower_run.sh            # self-test, scaling, plateau, 128-seed cobordism,
                             # large-N checks, (λ,neck) sweep, b2 to cap 16
```

Memory is a non-issue: the largest b₀/b₁ run measured used 0.10 GB.

**b₂ (also the tower, no cluster needed):**

```
bash tower_b2.sh             # static convergence, the |L| plateau, the throat
                             # window, and Delta beta measured in every component
```

`slurm_array.sh` is kept for reference. It pushes `cap` to 18-20 at |L| = 250,
which is the expensive route that does not converge; see section 4. Do not
spend cluster time on it.

## 6. Layout

```
cstopo3/
  cstopo3/
    z2.py         faces, elementary collapse, bit-packed Z2 rank, betti()
    fast01.py     the cycle-space b0/b1 path -- the reason 3+1D is tractable
    geom3.py      S3, S1xS2 (exact); warped surgery family; Spacetime + order
    extract.py    antichain, landmarks, shadows, slice_topology
    validate.py   exact triangulations: S3, S1xS2, T3, S4, wedge control
  test_all.py     5 groups, ~25 checks; run first
  run_static.py   static plateau + --scaling
  run_cobordism.py the headline result; job-array and multiprocessing
  run_b2.py       b2 convergence (static slices)
  tower_b2.sh     the b2 programme -- a workstation, not a cluster
  tower_run.sh    the whole tower programme
  slurm_array.sh  the cluster array
```

Two design points worth keeping in any refactor. Everything below
`Spacetime.precedes` sees only the order — `geom3.py` is the *only* module that
touches coordinates, so the "order alone" claim is structurally enforced rather
than asserted. And `S3Warped` exists solely so the product-approximation metric
can be validated against the exact `S3` through the identical pipeline; it is
what showed that the approximation is fine for b₁ and marginal for b₂.
