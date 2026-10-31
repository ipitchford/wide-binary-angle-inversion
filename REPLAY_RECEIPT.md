# Replay receipt

This is a producer-side replay, not an independent reproduction. It was run on 31 October 2026 on macOS (Apple silicon) with Python 3.11.14, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8, from a fresh copy of the final code, with the package directory placed on `sys.path` and the interpreter in isolated mode (`python3 -I -B`).

| Command | Result |
|---|---|
| `validate_and_plot.py --n-per-e 400000 --out OUT` | PASS. All checks passed. `results/summary.json`, `results/validation.json` and the four figure PNGs are byte-identical to the shipped copies. Runtime about 270 s. Pooled chi-square 119.98 on 115 degrees of freedom (p = 0.36). |
| `python3 -O ... validate_and_plot.py --n-per-e 400000 --out OUT-O` | PASS. `summary.json` identical to the normal run apart from `elapsed_seconds`. |
| `python3 -O ... negative_controls.py --out OUT-NC` | PASS. All 6 corruptions caught; `results/negative_controls.json` is byte-identical to the shipped copy. |
| `quick_checks.py --out DIR` (added for CI) | PASS in about 30 s under `-O`. |
| `latexmk -pdf manuscript.tex` | 21 pages, after inserting the DOI. Every page was inspected, and a PDF TeX preflight found no raw TeX. |

Publication packaging changed no code under test: the SHA-256 digests of `wide_binary_angles.py`, `validate_and_plot.py`, `negative_controls.py` and `requirements.txt` equal those of the replayed copy.

**Continuous integration.** GitHub Actions (`.github/workflows/replay.yml`) repeats on Ubuntu with Python 3.11: the checksum check, `quick_checks.py` and `negative_controls.py` under `python -O`, and a clean-tree check. `.github/workflows/full-replay.yml` runs the full suite on request; it does not block the release.
