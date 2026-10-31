# Agent navigation

## Identity

- Release slug: `wide-binary-angle-inversion`. Version `1.0.0-candidate`, tag `v1.0.0-candidate`.
- DOI 10.5281/zenodo.23838665; concept DOI 10.5281/zenodo.23838664.
- Repository: https://github.com/ipitchford/wide-binary-angle-inversion
- Release page: https://evidencepress.org/releases/wide-binary-angle-inversion/
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md)).
- Manuscript: [PDF](manuscript.pdf) (21 pages) and [LaTeX source](manuscript.tex). Figures: [likelihood_validation.png](figures/likelihood_validation.png), [moments_and_sensitivity.png](figures/moments_and_sensitivity.png), [powerlaw_index_information.png](figures/powerlaw_index_information.png), [inverse_instability.png](figures/inverse_instability.png).
- Overview and claims: [README.md](README.md), [CLAIMS.json](CLAIMS.json), [ASSURANCE.md](ASSURANCE.md), [VERIFICATION.md](VERIFICATION.md), [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
- Sources and provenance: [PROVENANCE.md](PROVENANCE.md), [SOURCES.md](SOURCES.md), [CITATION_AUDIT.md](CITATION_AUDIT.md), [RESEARCH_GATES.json](RESEARCH_GATES.json), [novelty_report.md](novelty_report.md), [target.yaml](target.yaml), [proof_core.md](proof_core.md).
- Citation and licences: [CITATION.cff](CITATION.cff), [.zenodo.json](.zenodo.json), [LICENSES.md](LICENSES.md), [PUBLIC_DOMAIN.md](PUBLIC_DOMAIN.md).

## Scientific map

**Setting.** A population of bound Kepler orbits with eccentricity law F on [0,1). Conditional on e: mean anomaly uniform, orbital frame Haar-random on SO(3), independent of phase. The observable Gamma in [0,pi] is the angle between the sky projections of relative position and relative velocity. Write c = |cos gamma|.

**Results.**
- Likelihood p(gamma|e) = [acosh(1/e) - 1{c>e} acosh(c/e)] / (pi sin^2 gamma). Theorem 3.1, **proved**; checked by quadrature and Kepler simulation ([wide_binary_angles.py](wide_binary_angles.py) `angle_pdf`).
- The angle law determines F (atoms included) by an explicit Abel inverse. Theorem 4.1, **proved** (`abel_inverse_left_cdf`).
- Power-law laws f(e) = (alpha+1) e^alpha give p_alpha(gamma) = B((alpha+1)/2, 1/2) (1 - c^(alpha+1)) / (2 pi sin^2 gamma); uniform angles iff thermal. Section 5, **proved**.
- E sqrt(1-e^2) = (pi/2) E sin^3 Gamma; every real-order circularity moment has an explicit observable. Theorem 6.1, **proved** (`fractional_moment_observable`).
- E e^2 = 1/2 + E cos 2Gamma - (1/2) E cos 4Gamma. Corollary 7.1, **proved** (`eccentricity_square_observable`).
- At most one symmetric unbiased single-angle observable per target. Theorem 6.2, **proved**.
- Power-law index estimators (`powerlaw_alpha_from_sin3`, `powerlaw_alpha_from_e2`), efficiency at least 0.9595 on the grid. Section 7.2, **formulas proved, numbers by quadrature**.
- Thermal invariance under azimuth-blind selection; separation binning exact under a log-uniform semimajor-axis law. Propositions 8.1-8.2, **proved**; simulated (Table 2).
- Inverse weakly continuous, discontinuous in total variation. Propositions 9.1-9.2, **proved**.

**Do not read:**
- any result as an analysis of real data: all data are synthetic;
- the angle of one binary as determining its eccentricity;
- uniform mean anomaly as uniform true anomaly;
- the grid maximum-likelihood fit as a consistent estimator (not proved);
- internal AI-assisted review as peer review.

## Code and replay

[quick_checks.py](quick_checks.py) runs the deterministic checks and reduced-size simulations in about 30 s. [negative_controls.py](negative_controls.py) runs 6 controls in about a minute. [validate_and_plot.py](validate_and_plot.py) `--n-per-e 400000 --out OUT` runs the full suite (about 5 minutes) and writes [results/summary.json](results/summary.json), [results/validation.json](results/validation.json) and the figures. Pinned versions: [requirements.txt](requirements.txt); environment: [ENVIRONMENT.txt](ENVIRONMENT.txt).

Replay record: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). CI: [.github/workflows/replay.yml](.github/workflows/replay.yml) (blocking) and [.github/workflows/full-replay.yml](.github/workflows/full-replay.yml) (on request).

## Useful next research question

Fit a Gaia wide-binary sample with the exact likelihood plus an explicit vector proper-motion error model and selection model (paper equation 49), and compare with simulation-grid fits. Theoretical: prove the 1/n total-variation decay, or a minimax rate for estimating F.
