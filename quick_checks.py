# SPDX-License-Identifier: MIT
"""Fast blocking checks for continuous integration (about a minute).

Run: python quick_checks.py --out DIR
Runs the deterministic checks (kernel normalisation, moment identities,
power-law mixtures, end-to-end Abel inversion, kernel derivative), an
independent Kepler simulation at reduced sample size, and the moment-identity
and estimator extensions at reduced size. Every check raises RuntimeError on
failure, also under python -O. Figures and results go to DIR, never into the
package. The full suite is validate_and_plot.py --n-per-e 400000.
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate_and_plot as vp  # noqa: E402


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out).resolve()
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    vp.ROOT = out
    rng = np.random.default_rng(20261028)
    det = vp.deterministic_checks()
    print("deterministic checks passed: max normalisation error",
          det["max_normalisation_error"], "max moment error", det["max_moment_error"])
    vp.physical_validation(50000, rng)
    print("Kepler simulation checks passed (n = 50,000 per eccentricity)")
    ext = vp.extensions(rng, n=20000, replicates=20, n_rep=2000)
    print("moment-identity extensions passed: max identity error", ext["max_identity_error"])
    (out / "results" / "quick_checks.json").write_text(json.dumps(dict(
        max_normalisation_error=det["max_normalisation_error"],
        max_moment_error=det["max_moment_error"],
        extension_max_identity_error=ext["max_identity_error"]), indent=2) + "\n")
    print("all quick checks passed")


if __name__ == "__main__":
    main()
