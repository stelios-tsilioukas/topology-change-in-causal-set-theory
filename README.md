# Topology change in Causal Set Theory

Code accompanying the paper

> **Topology change and the Gauss–Bonnet term in causal set theory:
> the sign of the cosmological term and the axion decay constant**
> S. A. Tsilioukas (University of Thessaly / National Observatory of Athens)

The paper argues that the Gauss–Bonnet term of a causal set is linear in the
element count, that the element count *is* the four-volume, and that the
coefficient of a cosmological constant term can therefore be read off without
forming any functional derivative. Two things in that argument are numerical
and are what this repository contains: the demonstration that a change of
spatial topology is recoverable **from the causal order alone**, and the
collection of checks behind the analytic claims.

The headline result:

```
S³  ──(index-1 surgery)──>  S¹×S²        Δβ = (0,+1,+1,0)
```

recovered from the order alone in **8 of 8 sprinklings at N = 2×10⁴**, in about
8 seconds per sprinkling on one core. Coordinates are used once, to generate the
causal relation, and never again.

---

## Layout

| directory | what it is | hardware |
|---|---|---|
| [`cstopo3/`](cstopo3/) | the 3+1D demonstration: the wormhole-sector transition | laptop; workstation for b₂ |
| [`cstopo/`](cstopo/) | the 2+1D demonstration: the trousers, where b₁ is directly meaningful | laptop |
| [`verification/`](verification/) | the one-off checks behind individual claims in the paper | laptop |

`cstopo3` is the one to read first. `cstopo` is the earlier and simpler
two-dimensional-slice version, kept because it is the control case in which the
answer is known independently.

## Install

```bash
git clone https://github.com/stelios-tsilioukas/topology-change-in-causal-set-theory.git
cd topology-change-in-causal-set-theory
pip install -r requirements.txt
```

`cstopo3` and `cstopo` need **numpy only**. Parts of `verification/` also use
SciPy, GUDHI and SymPy; see [requirements.txt](requirements.txt). Python 3.9 or
later. `mpi4py` is optional and is needed only to distribute seeds with `--mpi`;
every script runs without it.

## Run

```bash
cd cstopo3
python3 test_all.py                                       # ~2 min, run this first

python3 run_cobordism.py --N 20000 --seeds 8 --nL 250 \
        --cap 14 --procs 4                                # the headline result
```

Expected output:

```
  seed    1  before=[1, 0] after=[1, 1]  delta=[0, 1]  DETECTED  8s
  ...
  detection rate: 8/8
```

The 2+1D case is a single command and finishes in seconds:

```bash
cd cstopo && python3 test_all.py && python3 run_trousers.py --N 8000 --seeds 8
```

---

## What the pipeline does

```
sprinkle → causal relation → inextendible antichain A → landmarks L ⊂ A
        → nerve from past-shadows σ(y) ∩ L → GF(2) homology → Betti plateau
```

Everything below `Spacetime.precedes` sees only the order. `geom3.py` is the
only module that touches coordinates, so the "order alone" claim is enforced by
the module structure rather than merely asserted.

Two design points carry the 3+1D case, and both are explained at length in
[`cstopo3/README.md`](cstopo3/README.md):

**The resolution is fixed by the dimension of the slice, not by the sampling
density.** A Dowker/nerve complex reproduces the topology only when the cover
balls have radius `r ≈ 2s`, `s` the vertex spacing, so the facet size is
`m = ω_d (r/s)^d` — about 7 in two dimensions and about 14 in three. Carrying
the 2+1D value `c ≈ 10` over to a three-dimensional slice under-resolves it and
returns b₀ between 2 and 10² with b₁ of order 10²: a plateau of discreteness
artefacts rather than a topology. This, and not a shortage of elements, is why
the transition initially appeared to need `N ≳ 5×10⁴`.

**A cheap b₁.** Over Z₂ the class of a triangle boundary in the cycle space of
the 1-skeleton, in the fundamental basis of a spanning forest, is just the set
of its non-tree edges — so every column of ∂₂ has at most 3 nonzeros. Python
bigints beat numpy as the bitset at these sizes; that change alone was 30×.

## Cost

Cobordism, both slices, single core:

| N | ℓ | \|A\| before/after | \|L\| | wall | RSS | detected |
|---|---|---|---|---|---|---|
| 20 000 | 0.278 | 884 / 1 476 | 250 | 5.1 s | 0.05 GB | yes |
| 40 000 | 0.234 | 1 504 / 2 455 | 300 | 11.2 s | 0.06 GB | yes |
| 80 000 | 0.197 | 2 575 / 4 185 | 350 | 25.1 s | 0.08 GB | yes |
| 160 000 | 0.165 | 4 320 / 6 983 | 400 | 62.7 s | 0.10 GB | yes |

`|A| ∝ N^0.743`, as the N^(3/4) expectation requires. Memory is a non-issue for
b₀ and b₁: the largest run measured used 0.10 GB.

## What is *not* established

Stated plainly, because it constrains what the paper claims.

**b₂ is converged** (resolved after the first release; see
[`cstopo3/README.md`](cstopo3/README.md) §4). The full Betti vector is measured
on both sides of the transition,

    before = (1,0,0)    after = (1,1,1)    Δβ = (0,+1,+1)

at `N = 2×10⁵, |L| = 150, cap = 12, λ = 0.80, neck = 1.5`, in 206 s per seed on
one core and 0.31 GB. Poincaré duality is no longer a premise of
Δβ = (0,+1,+1,0); it is a consistency check the measured vectors satisfy. The
route that failed was raising `cap` at fixed `|L|`; the knob that works is a
*coarser* landmark set, plus a throat radius in the window λ ≈ 0.70–0.85, since
a thin throat has its 2-cycle filled in by the cover. Run
[`cstopo3/tower_b2.sh`](cstopo3/tower_b2.sh) — it needs a workstation, not a
cluster.

**The observable never uses a cross-split relation.** The before- and
after-slice witness windows lie entirely on their own side of `t_split`. What is
measured is the pair of slice topologies, which is exactly Δβ, but it is not a
test of the causal structure *through* the transition.

**The handle must not be short and fat.** λ = 0.6 with no neck gives 0/4: the
1-cycle has a representative shorter than the cover radius and is filled in.
λ ≤ 0.45 with a neck is robust. This is geometry rather than a bug, and the
physically relevant regime — a throat small compared with the region containing
it — is the one that works.

---

## Which script supports which claim

Section numbers refer to the current draft.

### Sec. II — the Euler density of a discrete space

| file | claim |
|---|---|
| `verification/kappa.py` | the identity χ = Σ_v κ(v); κ against the Regge deficit |
| `verification/kappa2.py` | κ in the stability window (alpha-complex); the locality test |
| `verification/kappa3.py` | flat torus → κ = 0 pointwise; an edge flip → an atomic defect |
| `verification/kappa4d.py` | the 4D case: S⁴ (χ = 2), Kühnel flat T⁴ (κ = 0 pointwise) |
| `verification/orderchi.py` | the order complex fails: trivial with a maximum, divergent otherwise |
| `verification/fcnerve.py` | the nerve of inclusive future cones fails |
| `verification/intnerve.py` | a bounded Alexandrov-interval cover fails |

### Sec. III — topology change and its detection

| file | claim |
|---|---|
| `cstopo3/` | **the 3+1D result**: S³ → S¹×S², 8/8 at N = 2×10⁴ |
| `cstopo/run_trousers.py` | **the 2+1D result**: (1,1) → (2,2), 8/8 at N = 8000 |
| `verification/homology.py` | the GF(2) Betti/Euler engine, validated on disk, S², T² |
| `verification/alpha2.py` | widest-plateau Betti numbers (the stability window) |
| `verification/confirm.py` | persistence bars on a sampled torus |
| `verification/cobordism2.py`, `slicefix.py` | the Morse bridge; the 1+1D pair of pants, χ(W) = −1 |
| `verification/bound.py` | **β = (1,1,1,1) does not identify S¹×S²**: the explicit 7-vertex wedge counterexample |
| `verification/nh.py`, `nh_exact.py`, `nh_check.py` | the minimal handle, n_h = 6, exhaustively |
| `verification/delta_beta_search.py` | search for the transition in small causal sets (`--nmax 18 --trials 200000`) |

### Sec. IV — the Gauss–Bonnet term and the scale hierarchy

| file | claim |
|---|---|
| `verification/paper_numbers.py` | **the authoritative reproduction** of every table and number in the paper |
| `verification/scales.py` | the scale hierarchy and the Poisson fluctuation table (see the note below) |
| `verification/rfix.py` | the 4D ordering fraction r = 0.502 and the abundance exponents |

### Sec. V — resolution, species, and the handle action

| file | claim |
|---|---|
| `verification/flowtest.py` | **the key negative result**: decimation does *not* destroy a handle by a power law. Detection stays at 1.00 down to f = 0.07, where the n_h = 6 power law predicts 1.7×10⁻⁶ |
| `verification/flowtest2.py` | the threshold is set by the handle's physical size (vary the tube radius) |
| `verification/gs.py` | the Giddings–Strominger derivation: sign, coefficient 3π²/4, charge-quantised spectrum |
| `verification/alpha_dep.py` | the S² cancellation; α enters only logarithmically |
| `verification/pathint.py` | the action-derived size spectrum |
| `verification/nariai.py`, `nariai2.py` | **the Nariai exclusion**: the flux-independent trace equation Λ = ½(a⁻²+b⁻²) > 0, and the b₂ − 2b₁ charge classification |
| `verification/bd_handle.py`, `bd_diag.py` | the direct Benincasa–Dowker handle action, and why it is out of reach: boundary terms, N ≳ 5×10⁴ |
| `verification/bd_decimate.py` | **the binomial thinning theorem** (Appendix B): decimation is exactly Sorkin smearing |
| `verification/bd_full.py` | the end-to-end check of that identity on a 4D sprinkling |
| `verification/higgs_fix.py` | the factor of 3 in the de Brito et al. triviality bound: their Eq. (37) is inconsistent with Eqs. (22) and (36) |
| `verification/higgs.py` | the bound as literally written, for comparison |
| `verification/coleman.py` | the dilute-instanton-gas rate ∂χ/∂V = Γ |

### Sec. VI and appendices

| file | claim |
|---|---|
| `verification/earlyde.py` | a deterministic Ω_DE ∝ H² reading is excluded: ΔN_eff ≈ 17, w = 0 in matter domination |
| `verification/dual.py` | **Appendix A**: order reversal does not enforce cancellation of the mean |
| `verification/csg.py`, `csg3.py` | classical sequential growth, exploratory |
| `verification/skeleton_checks.py` | assorted consistency checks on the argument structure |

---

## Two notes on reproducibility

**`paper_numbers.py` is authoritative.** It uses
Λ_obs ℓ_P² = 2.8497×10⁻¹²², and reproduces the paper's tables exactly:
S₁ = 285.95, f = 3.276×10¹⁶ GeV, a₀ = 6.215 ℓ_P, ℓ_def = 0.166 mm, and the exact
identity ℓ_def / ρ_Λ^(−1/4) = (4π)^(1/4) = 1.8828. Some of the earlier scripts —
`gs.py`, `alpha_dep.py`, `scales.py` — were written against a rounded
Λ ℓ_P² = 10⁻¹²² and so differ in the last digit (S₁ = 287.0 rather than 286.0,
ℓ_def = 0.215 mm rather than 0.166 mm). They are kept as written; where they
disagree with `paper_numbers.py`, `paper_numbers.py` is correct.

**`cstopo/run_wormhole.py` is superseded and is kept deliberately.** It is the
first 3+1D attempt, which carried the two-dimensional resolution over unchanged
and did not converge — modal β is noise at N ≤ 3×10⁴. It is retained so the
failure can be reproduced, since diagnosing it as a dimension-counting error
rather than a shortage of elements is what produced `cstopo3`. A useful
diagnostic, worth remembering: the elementary collapse removed essentially
nothing (172 310 → 172 256 simplices). A dense sample of a 3-manifold should
collapse by orders of magnitude, and that alone says the complex is not
manifold-like.

## Status of the long-running jobs

`cstopo3/tower_run.sh` is the full programme for a 10-core workstation (roughly
one to two hours): self-test, cost scaling, the resolution plateau, a 128-seed
cobordism, large-N checks, the (λ, neck) sweep, and b₂ to cap 16.
`cstopo3/tower_b2.sh` is the b₂ programme: static convergence, the `|L|`
plateau, the throat-radius window, and Δβ measured in every component. It also
fits on the same workstation. `slurm_array.sh` is kept for reference only — it
pushes `cap` to 18–20 at `|L| = 250`, which is the expensive route that does not
converge. No part of this repository needs a cluster.

## Reproducibility

Everything here runs on one core of an ordinary machine; no part of the
repository needs a cluster. Last verified on a clean checkout with Python 3.11
and the pinned minimum versions:

| suite | result |
|---|---|
| `cstopo3/test_all.py` | all tests pass, ~2 min |
| `cstopo/test_all.py` | all tests pass, seconds |
| `verification/` (40 scripts) | all run to completion |

Six of the verification scripts take minutes rather than seconds, by nature
rather than by defect: `bd_diag`, `bd_full` and `bd_handle` sprinkle and count
intervals at O(N²); `flowtest` and `flowtest2` sweep a persistence computation;
`delta_beta_search` is an exhaustive search over small causal sets.

`verification/paper_numbers.py` reproduces every table and derived number in the
paper and is the quickest way to check an installation end to end.

## Citing

If you use this code, please cite the accompanying paper. Machine-readable
metadata is in [CITATION.cff](CITATION.cff); GitHub renders it under *Cite this
repository* in the sidebar.

## License

MIT — see [LICENSE](LICENSE).
