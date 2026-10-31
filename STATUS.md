# Status

Version 1.0.0-candidate (31 October 2026). DOI [10.5281/zenodo.23838665](https://doi.org/10.5281/zenodo.23838665) (concept DOI [10.5281/zenodo.23838664](https://doi.org/10.5281/zenodo.23838664)). Scholarly creator: Anonymous.

**This is an unrefereed research candidate.** No journal and no named human specialist has reviewed it. The only checking was internal and AI-assisted (see [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)). That is not peer review, independent reproduction or formal verification.

**Target.** The question in Appendix B of Hwang, Ting and Zakamska, *MNRAS* 512 (2022) 3383 (doi:10.1093/mnras/stac675): does the distribution of the projected separation-velocity angle determine the eccentricity distribution of a wide-binary population?

Section and theorem numbers below refer to `manuscript.pdf`.

* **Proved in writing** (ideal model: uniform mean anomaly, Haar-random orientation, error-free projections)
  * Theorem 3.1: closed-form likelihood of the angle for one eccentricity.
  * Theorem 4.1: mixtures are Abel transforms; the angle law determines every eccentricity law, including atoms, by an explicit inverse.
  * Section 5: closed angle laws for power-law eccentricity populations; the thermal converse.
  * Theorem 6.1: every circularity moment of real order is the mean of an explicit bounded angle observable (for s >= 0, or when the moment is finite).
  * Theorem 6.2: an eccentricity functional has at most one symmetric unbiased single-angle observable.
  * Corollary 7.1: E e^2 and E e^4 as finite Fourier sums of the angle.
  * Proposition 8.1: thermal populations keep a uniform angle law under selections that ignore sky position angles (conditioned on semimajor axis and total mass).
  * Proposition 8.2: binning by projected separation is exact under a log-uniform (Opik) semimajor-axis law and otherwise reweights by an explicit factor.
  * Propositions 9.1 and 9.2: the inverse is discontinuous in total variation but continuous in the weak topology.
* **Formulas proved, numbers by quadrature**: Fisher information and estimator efficiencies for power-law indices (Section 7.2, Table 1).
* **Numerical demonstration only**: the grid maximum-likelihood fit (consistency and rate not proved); the observed 1/n decay in the total-variation example.
* **Open**: real-catalogue analysis with a selection and astrometric-error model; a proof of the 1/n decay; minimax rates; a closed joint angle-speed law; semiparametric efficiency.

**What is not claimed.**
* Any result about a real catalogue. All data are synthetic.
* Point identification of one binary's eccentricity from one angle.
* That angle-based eccentricity inference is new: it goes back to Tokovinin (1998) and was developed numerically by Tokovinin and Kiyaeva (2016), Tokovinin (2020) and Hwang et al. (2022).
* Historical priority. The literature search was bounded (`novelty_report.md`).
