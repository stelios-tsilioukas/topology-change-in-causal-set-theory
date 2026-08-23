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

**b₂ is not converged.** This is the honest negative result, and it constrains
what can be claimed.

`b₂` needs the 3-simplices and no cycle-space shortcut exists. It converges with
cap, but slowly and expensively (S³, N = 20 000, |L| = 250; correct value 0):

| cap | 10 | 12 | 14 | 16 |
|---|---|---|---|---|
| b₂ | 42 | 26 | 6 | 3 |
| raw 3-simplices | 68 270 | 171 395 | 383 321 | 710 921 |
| wall | 5 s | 19 s | 67 s | 198 s |

Clean trend to zero; cap ≈ 18–20 should reach it. **That is the cluster job.**

Counter-intuitively, a *denser* landmark set makes b₂ worse, because at fixed cap
more landmarks means a smaller `r/s` and a less solid complex (S³, cap 12):
`|L| = 150 → b₂ = 0` (correct), `250 → 15`, `400 → 54`.

At that sweet spot the **full** Betti vector comes out right, with exact metrics,
N = 40 000, |L| = 150, cap = 12: **S³ → (1,0,0) 4/4** and
**S¹×S² → (1,1,1) 4/4** (~55 s per slice).

**But not for the surgery model.** With the warped/approximate metric the same
parameters give b₂ ∈ {0,1,2} — unstable. So:

- **b₀, b₁: established.** Exact metrics, warped metrics, and the cobordism, all 8/8.
- **b₂: established for the exact static slices only.** Not for the surgery model
  at reachable resolution.

Since `b₂ = b₁` on a closed orientable 3-manifold by Poincaré duality, and
`b₃ = b₀`, **b₁ is already the discriminating observable** — the paper's
Δβ = (0,+1,+1,0) follows from Δb₁ = +1 plus duality. The paper reports it that way,
with the b₂ convergence as a separate, partially-completed check. A directly
measured (0,+1,+1,0) is not yet supported.

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

**On the cluster (32-cpu nodes) — b₂ only:**

```
sbatch slurm_array.sh        # 2 geometries x caps 14,16,18,20 x 8 seeds
python3 run_b2.py --merge results_b2/
```

Budget ~1 h wall on 32 tasks; the array is sized for 2 runs per task and cap 20
is the expensive one at ~1400 s. **Do not spend cluster time on b₀/b₁** — the
tower does 128 seeds in about two minutes.

If b₂ flattens above zero rather than reaching 0 and 1, the next knob is a
*smaller* |L| at the same cap, not a larger N.

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
  run_b2.py       b2 convergence; the cluster job
  tower_run.sh    the whole tower programme
  slurm_array.sh  the cluster array
```

Two design points worth keeping in any refactor. Everything below
`Spacetime.precedes` sees only the order — `geom3.py` is the *only* module that
touches coordinates, so the "order alone" claim is structurally enforced rather
than asserted. And `S3Warped` exists solely so the product-approximation metric
can be validated against the exact `S3` through the identical pipeline; it is
what showed that the approximation is fine for b₁ and marginal for b₂.
