# Assurance

This file sets out what has been checked and what has not. Each row is a separate kind of assurance, and no row implies another.

| Dimension | State | Evidence |
|---|---|---|
| Availability | Public GitHub release and Zenodo DOI 10.5281/zenodo.23838665 | Release assets, [SHA256SUMS](SHA256SUMS) |
| Written proofs | Complete proofs in the manuscript | `manuscript.pdf`, Sections 3-9 and Appendix A. Original proof notes for Theorems 3.1 and 4.1: `proof_core.md`. |
| Numerical checks | Passed; floating-point tolerances, not interval certificates | `validate_and_plot.py`: quadrature normalisation (max error 3.5e-14), moment identities (2.7e-14), end-to-end Abel inversion (3.1e-4 away from atoms), 5 x 400,000 Kepler simulations with chi-square and Kolmogorov-Smirnov tests at a family-wise level of 0.001, phase controls, selection simulations. See [VERIFICATION.md](VERIFICATION.md) and `results/`. |
| Optimised mode | Supported and replayed | Every check is an explicit test raising `RuntimeError`; the normal and `python -O` replays gave identical `summary.json` apart from elapsed time. |
| Semantic negative controls | Passed (6/6 caught) | `negative_controls.py`: kernel off by 1%, fractional observable scaled by 1.002, e^2 observable coefficient perturbed, inclination uniform in angle, uniform true anomaly, Abel inverse constant off by 5%. `results/negative_controls.json`. |
| Producer replay | Passed | Clean-copy replay of the final code in normal and `-O` modes. See [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). |
| Continuous integration | Fast blocking checks | `.github/workflows/replay.yml`: checksums, `quick_checks.py`, `negative_controls.py` under `python -O`, clean tree. The full 400,000-per-eccentricity suite runs separately and does not block. |
| Internal AI-assisted review | Done | Five-role review, an adversarial correctness review by a different model, a proof audit of new and rewritten proofs and a re-review. See [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). This is not peer review. |
| Implementation diversity | Partial | Reviewers recomputed the kernel, inverse, moment identities and instability numbers with their own code, including independent Kepler simulators. This is diversity within one producer workflow, not independent reproduction. |
| Independent reproduction by an unaffiliated party | Not established | |
| Specialist (human) review | Not established | |
| Editorial peer review | Not established | |
| Formal verification | Not attempted | No proof assistant was used. |
| Source fidelity | Checked | Hwang et al. (2022) Sections 2.1-2.5 and Appendices A-B read in full (arXiv:2111.01789v2); quotations checked. Tokovinin (1998) complete ADS scan inspected. See [SOURCES.md](SOURCES.md). |
| Novelty | Bounded search only | [RESEARCH_GATES.json](RESEARCH_GATES.json), `novelty_report.md` (28 October 2026, addendum 31 October 2026) |
| Citations | Every DOI resolved and its title matched | [CITATION_AUDIT.md](CITATION_AUDIT.md) |

## Scope and limits

* Everything assumes ideal sampling: uniform mean anomaly (or any phase law satisfying the pairing condition of Section 2), Haar-random orientation and error-free projections. Uniform mean anomaly is not uniform true anomaly.
* Selection by projected separation, proper-motion significance or quality cuts can change the phase and orientation laws; Section 8 identifies the exact cases that are safe.
* Measurement error needs its own model (equation 49 of the paper). Chance alignments with uniform angles are indistinguishable from a thermal component.
* Fine structure of the eccentricity law is nearly invisible in the angle law (Section 9).
* A replay shows that the shipped results are what the shipped programs produce. It does not by itself show that the programs implement the manuscript's definitions; that rests on the written proofs and the negative controls.
