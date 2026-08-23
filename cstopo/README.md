# cstopo — the 2+1D demonstration

> **Note.** `run_trousers.py` is live: it is the 2+1D result reported in the
> paper, where the spatial slice is a surface and b₁ is directly meaningful.
> `run_wormhole.py` is **superseded** by [`../cstopo3/`](../cstopo3/) and is kept
> only so that its failure can be reproduced — it carried the two-dimensional
> resolution over to a three-dimensional slice unchanged, which under-resolves
> it. See [`../cstopo3/README.md`](../cstopo3/README.md) §1.


Dependency-light (numpy only), job-array parallel, MPI hook for large runs.

## Install / run

```bash
pip install numpy          # that is all
python run_trousers.py --N 8000 --seeds 8
```

## Two demonstrations

### 1. `run_trousers.py` — 2+1D, **works on a laptop**

Spatial slices go one circle → two circles: the 2+1D analogue of
S³ → S¹×S² with b₁ directly meaningful. Verified:

| N | wall time | early slices | late slices | detection |
|---|---|---|---|---|
| 6 000 | 1.7 s | (1,1) | (2,2) | — |
| 8 000 | 3.8 s | (1,1) | (2,2) | **8/8 seeds** |
| 12 000 | 4.9 s | (1,1) | (2,2) | — |

Plateau stable over shadow-caps 6–12. **No cluster needed**; a cluster only
buys ensemble statistics, which is embarrassingly parallel.

### 2. `run_wormhole.py` — 3+1D Giddings–Strominger, **needs a cluster**

Target: Δβ = (0,+1,+1,0), hence δχ = −2.

**Run `--scaling` first.** Current status at accessible N:

```
       N    |A|   plateau width       modal beta     time
    3000     44             4.0    (1, 11, 0, 0)     0.9s
   10000    172             4.0    (1,  0,29, 1)     6.4s
   30000    509             4.0   (1, 10,36, 0)    48.2s
```

Not converging: modal β is noise, and the cost is ~O(N²) (10× in N → 50× in
time). Extrapolating, N = 10⁵ is ~10 min/seed and N = 10⁶ is ~17 h/seed with a
~125 GB causal matrix. **Do not request a large allocation before the scaling
study shows the plateau stabilising.** A negative result at accessible N is a
statement about the observable, not about the physics — report the scaling.

## Pipeline (order-only)

```
sprinkle → causal matrix → inextendible antichain → MRS nerve
        → GF(2) homology → Betti plateau → Morse count Σ(−1)^λ
```

After the causal matrix is built the embedding is discarded: everything
downstream uses only the order, which is the point of the demonstration.

## Design notes

**Thickening is by shadow cardinality, not layer count.** An element far to the
future of an antichain has essentially all of it in its past, so its shadow
carries no local information and the nerve degenerates. The resolution
parameter is `shadow_cap`, the maximum cardinality of a cover set — the
discrete analogue of a small ball. This was the single fix that made the
demonstration work.

**Plateau, not widest plateau.** The correct Betti vector need not occupy the
widest plateau; a competing plateau exists where components have merged. The
scale must be chosen physically (ℓ_cg), not by maximising width. Scan `--caps`
and inspect.

**β = (1,1,1,1) is necessary but not sufficient** for S¹×S²: the wedge
S¹∨S²∨S³ has the same Betti vector. Use `--pseudomanifold` for the cheap
necessary check (every 2-face in exactly two 3-simplices).

## Cluster usage

```bash
#SBATCH --array=0-31
python run_trousers.py --N 40000 --seeds 128 \
    --task $SLURM_ARRAY_TASK_ID --ntasks 32 --out results/
python run_trousers.py --merge results/
```

MPI alternative: `mpirun -n 32 python run_trousers.py --mpi ...` (needs mpi4py).

## Files

- `cstopo/causet.py` — packed-bit causal matrix (N²/8 bytes), sprinkling,
  light-cone-windowed relation building. Verified against brute force.
- `cstopo/geometry.py` — trousers, cylinder (control), GS wormhole.
- `cstopo/topology.py` — antichains, layers, MRS nerve, GF(2) homology
  (validated on point/interval/circle/disk/S²/T²/two-circles), stability scan,
  Morse count.
