# Figure data

Output of `figdata.py`, and the input from which Fig. 2 of the paper is drawn.

| file | contents |
|---|---|
| `fig_plateau.json` | (b₀, b₁), landmark count and witness count against the shadow-cardinality resolution `c`, for S³ and S¹×S², 8 seeds at N = 2×10⁴, \|L\| = 250. Panels (a) and (b). |
| `fig_b2.json` | the b₂ convergence in `c` for both geometries. Panel (c), the unfinished part. |

Regenerate with `python3 figdata.py` from the `cstopo3/` directory. It writes to
the working directory, so copy the results here to replace these.
