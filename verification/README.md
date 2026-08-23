# Verification scripts

One-off checks supporting individual claims in the paper. Each is standalone
and prints its own result; run any of them with `python3 <name>.py`. The
directory is deliberately flat so that the cross-imports between scripts
(`homology.py`, `fcnerve.py`, `nh.py`, `flowtest.py`, `bd_decimate.py`,
`kappa2.py` are imported by others) keep working.

The mapping from script to paper claim is in the
[top-level README](../README.md#which-script-supports-which-claim).

## Start here

    python3 paper_numbers.py     # reproduces every table and number in the paper

## Notes

`paper_numbers.py` is authoritative. `gs.py`, `alpha_dep.py` and `scales.py`
were written against a rounded Λ ℓ_P² = 10⁻¹²² and differ from the paper in the
last digit; they are kept as written.

Four scripts are slow (minutes, not seconds), by nature rather than by defect:

| script | why |
|---|---|
| `bd_diag.py`, `bd_full.py` | sprinkle and count intervals, O(N²) |
| `flowtest2.py` | a tube-radius sweep, each point a full persistence computation |
| `delta_beta_search.py` | an exhaustive search; takes `--nmax` and `--trials` |

`handle_validation.py` is not included: it depended on `ripser`, which no longer
builds on current Python. Its function is covered by `cstopo3/test_all.py`.

## Dependencies

numpy throughout; scipy and GUDHI for some. See
[requirements.txt](../requirements.txt).
