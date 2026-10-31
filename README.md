# An exact projected-angle likelihood and population inversion for wide-binary eccentricities

Evidence Press research candidate, version 1.0.0-candidate, 31 October 2026.
Author: Anonymous. Text CC0-1.0; code MIT.

- DOI: [10.5281/zenodo.23838665](https://doi.org/10.5281/zenodo.23838665) (all versions: [10.5281/zenodo.23838664](https://doi.org/10.5281/zenodo.23838664))
- Release page: https://evidencepress.org/releases/wide-binary-angle-inversion/
- Status: unrefereed candidate; internal AI-assisted review only ([STATUS.md](STATUS.md), [ASSURANCE.md](ASSURANCE.md))
- For agents: [AI_INDEX.md](AI_INDEX.md) and [CLAIMS.json](CLAIMS.json)

## What this is about

A wide binary star takes thousands of years to orbit, so its eccentricity
usually cannot be found by following the orbit. Astronomers instead measure,
for each pair, the angle `Gamma` between the projected separation and the
projected relative velocity, and study the distribution of that angle over
many pairs. This package answers four questions about that method under the
ideal model (random orbital phase, uniformly random 3D orientation):

1. **What is the exact likelihood of one angle?** A closed formula
   (Theorem 3.1): `p(gamma|e) = [acosh(1/e) - 1{|cos gamma|>e} acosh(|cos gamma|/e)] / (pi sin^2 gamma)`.
2. **Does the angle distribution determine the eccentricity distribution?**
   Yes, for every eccentricity law, including atoms, by an explicit Abel
   inverse (Theorem 4.1). This answers the uniqueness question raised in
   Appendix B of Hwang, Ting and Zakamska (2022). The inverse is continuous in
   the weak topology (Proposition 9.2) but not in total variation
   (Proposition 9.1), so fine detail is nearly invisible.
3. **What can be estimated directly, without inversion?** Every moment
   `E (1-e^2)^(s/2)` is the mean of an explicit observable (Theorem 6.1), e.g.
   `E sqrt(1-e^2) = (pi/2) E sin^3 Gamma` and
   `E e^2 = 1/2 + E cos 2Gamma - (1/2) E cos 4Gamma` (Corollary 7.1). These are
   the only symmetric unbiased single-angle observables for their targets
   (Theorem 6.2). Power-law indices follow in one line, with near-optimal
   efficiency (Section 7.2, Table 1).
4. **Which observational selections keep the model exact?** Populations that
   are thermal at each semimajor axis keep a uniform angle law under any
   selection on true lengths and line-of-sight angles that ignores sky position
   angles (Proposition 8.1). Binning by projected separation is exact under a
   scale-free log-uniform semimajor-axis law and otherwise reweights by an
   explicit factor (Proposition 8.2).

No observational catalogue is analysed. Real-catalogue use needs a selection
and astrometric-error model (equation 49 of the paper).

## Files

- `manuscript.pdf` / `manuscript.tex`: the paper, with complete proofs.
- `wide_binary_angles.py`: likelihoods, observables, estimators, simulators.
- `validate_and_plot.py`: every numerical check, table and figure.
- `negative_controls.py`: shows that corrupted kernels, observables,
  simulators and inverse constants are caught.
- `results/`: `summary.json`, `validation.json`, `negative_controls.json`.
- `figures/`: four figures (PDF and PNG).
- `proof_core.md`: the original proof notes for Theorems 3.1 and 4.1.
- `novelty_report.md`: the bounded prior-art search, with a 31 October addendum.
- `VERIFICATION.md`: what was checked and what was not.
- `REVIEW_SUMMARY.md`: the internal AI-assisted review, its objections and how
  each was handled.
- `quick_checks.py`: the fast checks run by continuous integration (about 30 s).
- Package records: `STATUS.md`, `ASSURANCE.md`, `PROVENANCE.md`, `REPLAY_RECEIPT.md`,
  `SOURCES.md`, `CITATION_AUDIT.md`, `RESEARCH_GATES.json`, `ENVIRONMENT.txt`,
  `CITATION.cff`, `.zenodo.json`, `LICENSES.md`, `SHA256SUMS`.

## Reproduce

Python 3 with NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8:

```bash
python validate_and_plot.py --n-per-e 400000 --out replay
python negative_controls.py --out replay-nc
```

The seed is 20261028. The first command writes `replay/results` and
`replay/figures` and leaves the shipped results unchanged; it takes a few
minutes. Every check raises an error on failure, including under `python -O`.
All data are synthetic.

To build the paper: `latexmk -pdf manuscript.tex`.

## Use

```python
import numpy as np
from wide_binary_angles import (angle_pdf, powerlaw_angle_pdf, moment_observable,
    eccentricity_square_observable, powerlaw_alpha_from_sin3)

angle_pdf(np.radians([30, 60, 90]), 0.6)        # likelihood per radian on [0, pi]
# For angles from an ideal sample:
# np.mean(moment_observable(angles, 1))          # estimates E sqrt(1-e^2)
# np.mean(eccentricity_square_observable(angles)) # estimates E e^2
# powerlaw_alpha_from_sin3(angles)              # alpha for f(e) = (alpha+1) e^alpha
```

## Scope that must be kept

Uniform mean anomaly is not uniform true anomaly. Selection by projected
separation, proper-motion significance or quality cuts can change the phase
and orientation laws (Section 8). Measurement error needs its own model.
Chance alignments with uniform angles are indistinguishable from a thermal
component. The grid maximum-likelihood fit is a demonstration; its consistency
and rate are not proved. The work is unrefereed; checking was internal and
AI-assisted.

## Cite

Anonymous (2026). *An exact projected-angle likelihood and population inversion
for wide-binary eccentricities* (version 1.0.0-candidate). Evidence Press.
https://doi.org/10.5281/zenodo.23838665
