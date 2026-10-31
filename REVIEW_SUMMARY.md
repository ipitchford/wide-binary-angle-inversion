# Review summary (public, rating-free)

This candidate received an internal, AI-assisted review on 31 October 2026.
It is not external peer review. Reports were produced independently in five
roles (editor, methodology, astronomy domain, cross-disciplinary perspective,
devil's advocate), plus one adversarial correctness review, a proof audit of
the new and rewritten proofs, and a re-review of the revision. Routes: Claude
Code subagents (devil's advocate, methodology, domain); OpenAI gpt-5.6-sol
(adversarial correctness review of the submitted version); Codex gpt-6-astra
(perspective, editor, proof audit, re-review). Different AI reviewers agreeing
is not independent proof.

## Correctness

No reviewer found an error in the main theorems. Several recomputed the
kernel, the inverse, the moment identities and the instability numbers with
their own code, including independent Kepler simulators. The adversarial
review asked for a formal treatment of two steps: the independence of the
orbital frame and the phase, and the choice of a canonical representative of
the angle density before differentiating. Both are now written out (Sections
3.1 and 4). The proof audit of the revision found seven repairable gaps of
scope (negative-order moments, the class in which unbiased observables are
unique, a moment-existence condition for Makarov's family, the contamination
estimator, the conditioning variables in the thermal selection result, an
improper semimajor-axis measure, and wording of statistical claims). All were
repaired. The re-review found a further error in the numerical demonstration:
the grid fit's distribution function had been interpolated rather than
evaluated as a step function. It was corrected and the numbers regenerated.

## Significance and framing

- Hwang et al.'s Appendix B raises uniqueness as "a related question" and
  suspects approximate uniqueness. The paper now quotes it and answers both
  parts: exact uniqueness holds; approximate uniqueness holds weakly
  (Proposition 9.2, added) and fails in total variation.
- "Nonparametric inversion" suggested an estimator that did not exist. The
  title now says "population inversion", and an approximate grid
  maximum-likelihood fit is demonstrated (consistency not proved).
- Reviewers asked what practitioners can use. Added: estimators of the mean
  square eccentricity and of power-law indices in both common
  parametrisations; Fisher information per angle and estimator efficiencies;
  the general dual-observable theorem.
- The astronomy reviewer and devil's advocate showed that the original
  selection experiment used an unobservable cut. Added: an observation model;
  an exact thermal-invariance result for selections that ignore sky position
  angles, which proves a numerical finding of Hwang et al.; an exact
  separation-binning result (invariance under a scale-free log-uniform
  semimajor-axis law, explicit reweighting otherwise) with simulations for
  realistic bins; a note that uniform contaminants mimic a thermal component.
- Missing literature was added from Crossref records: Benisty et al. (2022),
  Jeans (1919), Shatsky (2001), Hamilton (2022), Makarov parts II and III,
  El-Badry et al. (2021), Hwang et al. (2022, twins), Wicksell (1925),
  Gorenflo and Vessella (1991), Groeneboom and Jongbloed (1995), Hausdorff
  (1923), Bhatia and Davis (2000).

## Validation

The methodology reviewer showed that one advertised check (atomic inversion)
could not fail, that every check was an `assert` and vanished under
`python -O`, and that the distributional tests did not gate the run. All three
are fixed: an end-to-end inversion from the angle density, explicit checks,
Kolmogorov-Smirnov tests with a family-wise gate, phase controls (uniform
eccentric anomaly must pass, uniform true anomaly must fail), an `--out`
option, and a negative-control script.

## Not done (open)

- No real catalogue was analysed; that needs a vector proper-motion error
  model and a selection model.
- No proof of the observed 1/n decay of the total-variation example, no
  minimax rate, and no consistency theory for the grid fit.
- No closed joint law for angle and speed.
- No semiparametric efficiency theorem.
- The "rigid-angle projection" literatures (jet-disc and spin-velocity
  alignments) were not searched for a second route to uniqueness.
- Hoeffding intervals remain conservative; sharper intervals are future work.
