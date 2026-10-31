"""Semantic negative controls: deliberately corrupted ingredients must be caught.

Run: python negative_controls.py --out DIR
Each control patches one ingredient inside validate_and_plot's namespace and
runs the relevant check. The control passes only if the check raises
RuntimeError. Works under python -O (no assert statements are used).
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import validate_and_plot as vp  # noqa: E402
import wide_binary_angles as wba  # noqa: E402


def wrong_kernel(gamma, eccentricity):
    # Plausible algebra slip: acosh(c/e) replaced by acosh(c/e) at 1% larger e.
    e = np.minimum(np.asarray(eccentricity, float) * 1.01, .99999)
    return wba.angle_pdf(gamma, e)


def wrong_moment(gamma, s):
    return wba.fractional_moment_observable(gamma, s) * (1 + 0.002)


def wrong_e2(gamma):
    x = np.cos(2 * np.asarray(gamma, float))
    return 1 + x - 0.98 * x * x


def uniform_inclination_sampler(eccentricity, rng, phase="mean", return_state=False):
    # Orientation error: inclination uniform in angle instead of uniform cos(i).
    e = np.asarray(eccentricity, float)
    m = rng.uniform(0, 2*np.pi, e.shape)
    anom = m.copy()
    for _ in range(80):
        anom = anom - (anom-e*np.sin(anom)-m)/(1-e*np.cos(anom))
    j = np.sqrt(1-e*e)
    rx, ry = np.cos(anom)-e, j*np.sin(anom)
    vx, vy = -np.sin(anom), j*np.cos(anom)
    om = rng.uniform(0, 2*np.pi, e.shape)
    ci = np.cos(rng.uniform(0, np.pi, e.shape))
    cw, sw = np.cos(om), np.sin(om)
    r1, r2 = rx*cw-ry*sw, ci*(rx*sw+ry*cw)
    v1, v2 = vx*cw-vy*sw, ci*(vx*sw+vy*cw)
    g = np.arctan2(np.abs(r1*v2-r2*v1), r1*v1+r2*v2)
    if return_state:
        return dict(gamma=g, cos_inclination=ci, projected_separation=np.hypot(r1, r2),
                    kepler_residual=0.0)
    return g


def true_anomaly_sampler(eccentricity, rng, phase="mean", return_state=False):
    return wba.sample_kepler(eccentricity, rng, phase="true", return_state=return_state)


def wrong_abel(angle_density, r, singular_points=(), h=1e-6):
    return 1.05 * wba.abel_inverse_left_cdf(angle_density, r, singular_points, h)


CONTROLS = [
    ("kernel off by 1% in eccentricity", "angle_pdf", wrong_kernel,
     lambda rng: vp.deterministic_checks()),
    ("fractional observable scaled by 1.002", "fractional_moment_observable", wrong_moment,
     lambda rng: vp.extensions(rng, n=20000, replicates=20, n_rep=2000)),
    ("e^2 observable coefficient perturbed", "eccentricity_square_observable", wrong_e2,
     lambda rng: vp.extensions(rng, n=20000, replicates=20, n_rep=2000)),
    ("inclination uniform in angle", "sample_kepler", uniform_inclination_sampler,
     lambda rng: vp.physical_validation(100000, rng)),
    ("uniform true anomaly", "sample_kepler", true_anomaly_sampler,
     lambda rng: vp.physical_validation(100000, rng)),
    ("Abel inverse constant off by 5%", "abel_inverse_left_cdf", wrong_abel,
     lambda rng: vp.deterministic_checks()),
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    out = Path(args.out).resolve()
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "results").mkdir(parents=True, exist_ok=True)
    vp.ROOT = out
    rows = []
    for label, name, replacement, run in CONTROLS:
        original = getattr(vp, name)
        setattr(vp, name, replacement)
        try:
            run(np.random.default_rng(7))
            caught, message = False, "no error raised"
        except RuntimeError as exc:
            caught, message = True, str(exc)
        finally:
            setattr(vp, name, original)
        rows.append(dict(control=label, patched=name, caught=caught, message=message))
        print(f"{'CAUGHT' if caught else 'MISSED'}: {label} -> {message}")
    (out / "results" / "negative_controls.json").write_text(json.dumps(rows, indent=2) + "\n")
    if not all(r["caught"] for r in rows):
        raise SystemExit("negative control missed")
    print("all negative controls caught")


if __name__ == "__main__":
    main()
