# Provenance

- **Scholarly creator.** Anonymous, under the standing Evidence Press protocol. Ian Pitchford maintains and publishes the repository; that role is not authorship.
- **Production.** The research, proofs, code, manuscript and revisions were produced with AI language-model tools under human direction, as the manuscript's status statement says. [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md) records the review history.
- **Timeline (2026).**
  - 28 October: target definition (`target.yaml`), bounded novelty search (`novelty_report.md`), proof notes (`proof_core.md`), first manuscript and code.
  - 31 October: internal five-role review, adversarial correctness review and clean replay; one revision round adding the selection results, dual observables, Fourier eccentricity moments, power-law estimators with Fisher information and weak continuity; proof audit, re-review and final fixes; clean replay of the final code in normal and `-O` modes with negative controls; publication packaging.
- **Publication-only changes.** The version, DOI and repository URL were added to the title block and the reproducibility appendix, and the PDF was rebuilt. Package records (this file, `STATUS.md`, `ASSURANCE.md`, `AI_INDEX.md`, `CLAIMS.json`, citation and licence files), a fast CI script (`quick_checks.py`) and CI workflows were added, and the README gained an identity header and pointers to them. No scientific content, code under test or result changed in packaging.
- **Not shipped.** The private review reports and the editor's working notes are kept privately. Their substantive objections and resolutions are in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md).
