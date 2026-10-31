# Verification record

Date: 31 October 2026. Candidate: "An exact projected-angle likelihood and
population inversion for wide-binary eccentricities".

## What is established, and how

| Claim | Basis | Numerical support |
|---|---|---|
| Exact likelihood (Thm 3.1) | Written proof | Quadrature normalisation at 8 eccentricities (max discrepancy 3.5e-14); 5 x 400,000 independent Kepler simulations, chi-square and Kolmogorov-Smirnov tests; phase controls |
| Measure-level inverse (Thm 4.1) | Written proof | End-to-end inversion from the angle density for an atomic law and three power laws (max error 3.1e-4 away from atoms) |
| Closed power-law laws (Sec. 5) | Written proof | Mixture quadrature at 12 points |
| Real-order moments (Thm 6.1), bounds | Written proof | Identities at 6 orders x 5 eccentricities (max 2.7e-14); dense-grid range check; Kepler simulation |
| Dual observables (Thm 6.2) | Written proof | Worked examples only |
| Even eccentricity moments (Cor 7.1) | Written proof | Quadrature and Kepler simulation |
| Power-law index estimators, Fisher information | Formulas proved; information and efficiencies by quadrature | Table 1; 501-point efficiency grid (min 0.9595 at alpha = 0.5); 400 replicate thermal samples |
| Thermal selection invariance (Prop 8.1) | Written proof | 2,000,000-binary simulation with separation bin and speed cut |
| Separation binning (Prop 8.2) | Written proof | Six simulated populations (Table 2) |
| TV discontinuity (Prop 9.1) | Written proof | Fresnel evaluation, n = 1..64 |
| Weak continuity (Prop 9.2) | Written proof | none |
| Grid maximum-likelihood fit | Demonstration only | Step-CDF errors 0.103, 0.061, 0.036 at N = 2,000, 20,000, 100,000 |

Numerical tolerances are floating-point checks, not interval certificates.
Proofs are written arguments; no proof assistant was used.

## Replay

Command (from a fresh copy, outside the package directory):

```bash
python validate_and_plot.py --n-per-e 400000 --out OUT
python negative_controls.py --out OUT-NC
```

Environment of the recorded run: Python 3.11.14, NumPy 2.3.5, SciPy 1.17.0,
Matplotlib 3.10.8 (macOS arm64). Seed 20261028. Runtime about 270 s on a
loaded machine. The run was made in normal and optimised (`-O`) mode; the two
`summary.json` files are identical apart from elapsed time. Every check is an
explicit test that raises an error. Distributional tests are gated at a
family-wise level of 0.001. All six negative controls (wrong kernel, wrong
observables, wrong inclination law, uniform true anomaly, wrong Abel constant)
are caught; see `results/negative_controls.json`.

An earlier replay of the original 28 October package on the same platform
reproduced its recorded numbers to floating-point precision.

## Review

Internal AI-assisted review only; see `REVIEW_SUMMARY.md`. No external peer
review.

## Not established

- Any result about a real catalogue; no observational data are included.
- Validity of the ideal phase and orientation assumptions for any selected
  sample beyond Propositions 8.1 and 8.2.
- A rate for the 1/n attenuation, a minimax rate, or consistency of the grid fit.
- Semiparametric efficiency.
- Exhaustive prior-art exclusion (see `novelty_report.md`).
