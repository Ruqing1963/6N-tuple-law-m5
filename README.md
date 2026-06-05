# The Tuple-Size Law at m=5: Resolving the Sparse Regime (Part XIX)

Part XVIII confirmed `d ln rho_m / d ell -> -m` at m=4 and predicted the prime
quintuplet (m=5): slope **-5**, ratio **5:2** against the twin. This is that test
— and the rung where the sparse regime the law had been approaching finally has
to be met head-on. **The prediction holds, and the manner of its holding is the
story.**

**Sparsity bites (S9->S10).** With only 2,936 quintuplets on S9 and 16,570 on
S10, the slope reads **-5.43**, ratio **2.62** — overshooting the `2.5 +/- 0.06`
bound, exactly the caveat Part XVIII raised (large m=5 tail + sparsity noise both
pull the estimate up).

**Depth resolves it (S10->S11).** Recomputing from sieve primitives gives
**102,254** quintuplets on S11 (twin counted in the same pass: **196,963,369**),
so the deeper, better-sampled slope is **-5.065 +/- 0.09**, ratio **2.4488 +/-
0.042** — inside the bound, on the law. The tail shrank from +0.43 to +0.065
(deeper = smaller tail); the theory predicted the bias and the cure, and the
computation obeyed: `-5.43 -> -5.065`, `2.62 -> 2.449`.

| pattern | m | \|slope\| (deepest reliable) | ratio to twin (target) |
|---|---:|---:|---:|
| twin | 2 | 2.072 (S9->S10) | 1.0000 |
| triplet | 3 | 3.111 (S9->S10) | 1.5011 (3/2) |
| quadruplet | 4 | 4.175 (S9->S10) | 2.0145 (2) |
| **quintuplet** | **5** | **5.065 (S10->S11)** | **2.4488 (5/2)** |

> **Four-point law.** `|d ln rho_m / d ell| = kappa m`; ratios
> `1 : 1.5011 : 2.0145 : 2.4488` against `1 : 3/2 : 2 : 5/2`. Four points through
> the origin fix it: tuple size, and nothing else, sets the decay rate.
>
> **The sparse-regime instruments, both validated.** The count-aware estimator
> (Part XVIII) was stress-tested in the sparsest data yet and held; the residual
> bias it could not remove was removed by descent to S11. Weighting and depth are
> now both characterised by where each is needed.
>
> **Scope.** First moments only; the omega-variance / Erdos-Kac law on twin
> centres is open (Part XII). Hardy-Littlewood inputs; rho is the Part I scale. No
> infinitude is claimed. The quintuplet slope is at a deeper baseline (S10->S11)
> than the others (S9->S10); since the twin shifts only -2.072 -> -2.068 between
> them, the comparison is fair to <0.5%, and the 2.4488 ratio uses twin and
> quintuplet from the same S11 pass.

**This closes the empirical mining of the series.** m=6 (the prime sextuplet
`(6N+1, 6N+5, 6N+7, 6N+11, 6N+13, 6N+17)`) must give slope **-6** and ratio
**3:1** against the twin; it needs shells beyond S11 (~10x rarer again) and is
left as the open challenge for supercomputer centres. A twin ratio outside
`3 +/- 0.07` would falsify the law.

Part I: doi:10.5281/zenodo.20470367 · XVII: doi:10.5281/zenodo.20541575 ·
XVIII: doi:10.5281/zenodo.20546249

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── .zenodo.json
├── data/
│   ├── tuple_flow.csv          shell, ell, and (count, rho, ln_rho) for all four patterns (twin to S11)
│   └── quintuplet_counts.csv   shell, quint_count, centres, rho_quint, source (S2-S11)
├── code/
│   ├── quintuplet_geodynamics.py   four-pattern verdict + final figure (count-aware)
│   └── quintuplet_deep.py          same-pass twin+quintuplet from-scratch recompute (S10, S11)
├── figures/                fig_paper19_tuplelaw_m5.{pdf,png}
└── paper/                  Chen_6N_Paper19.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# Four-pattern verdict + the final figure (uses the embedded counts).
python code/quintuplet_geodynamics.py

# Same-pass twin + quintuplet recompute from scratch (the bedrock check).
python code/quintuplet_deep.py                 # S10 and S11
SHELLS=11 python code/quintuplet_deep.py        # S11 only (~1-2 hours)
```

On Windows PowerShell set env vars separately, e.g.
`$env:SHELLS="11"; python code/quintuplet_deep.py`.

### Conventions (same as Parts I-XVIII)

- Skeleton 6N±1; rho_m = #{m-tuple centres} / #{centres N} per shell;
  centres per shell = 1.5e(k-1); ell = ln ln X.
- quintuplet (6N-1, 6N+1, 6N+5, 6N+7, 6N+11) = pattern (0,2,6,8,12).
- Engine: segmented value sieve; every count recomputed from primitives.

## License

MIT — see `LICENSE`.
