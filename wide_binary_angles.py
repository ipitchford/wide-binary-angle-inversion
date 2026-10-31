"""Exact projected position--velocity angle laws for ideal Keplerian binaries.

Angles are in radians. Eccentricities satisfy 0 <= e < 1. The model assumes
uniform mean anomaly and an independent Haar-uniform orbital orientation.
These functions do not correct survey selection or astrometric error.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import beta, digamma


def angle_pdf(gamma, eccentricity):
    """Density on [0, pi], with stable evaluation at small projected angles.

    Arguments follow NumPy broadcasting. An exact circular orbit has an
    integrable infinite density at pi/2. Endpoint values are limiting values.
    """
    g, e = np.broadcast_arrays(np.asarray(gamma, float),
                               np.asarray(eccentricity, float))
    if np.any((e < 0) | (e >= 1) | ~np.isfinite(e)):
        raise ValueError("eccentricity must be finite and in [0, 1)")
    if np.any((g < 0) | (g > np.pi) | ~np.isfinite(g)):
        raise ValueError("gamma must be finite and in [0, pi]")
    s2 = np.sin(g) ** 2
    c = np.abs(np.cos(g))
    j = np.sqrt((1 - e) * (1 + e))
    out = np.empty(g.shape)
    ends = (g == 0) | (g == np.pi)
    out[ends] = 1 / (2 * np.pi * j[ends])
    singular = (e == 0) & (g == np.pi / 2)
    out[singular] = np.inf
    inner = (c < e) & ~ends & ~singular
    # acosh(1/e) = log((1 + sqrt(1-e^2))/e)
    out[inner] = (np.log1p(j[inner]) - np.log(e[inner])) / (
        np.pi * s2[inner])
    outer = ~(ends | singular | inner)
    co, eo, jo, so = c[outer], e[outer], j[outer], s2[outer]
    d = np.sqrt(np.maximum((co - eo) * (co + eo), 0))
    # Stable acosh difference: log((1+j)/(c+d)); rationalise j-d.
    z = so * (1 / (1 + co) + 1 / (jo + d)) / (co + d)
    out[outer] = np.log1p(z) / (np.pi * so)
    return float(out) if out.ndim == 0 else out


def powerlaw_angle_pdf(gamma, alpha):
    """Exact density for f(e)=(alpha+1)e**alpha, alpha > -1."""
    if not np.isfinite(alpha) or alpha <= -1:
        raise ValueError("alpha must be finite and greater than -1")
    g = np.asarray(gamma, float)
    if np.any((g < 0) | (g > np.pi) | ~np.isfinite(g)):
        raise ValueError("gamma must be finite and in [0, pi]")
    q = alpha + 1
    s2 = np.sin(g) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        logc = np.where(s2 < .5, .5*np.log1p(-s2), np.log(np.abs(np.cos(g))))
        logc = np.where(g == np.pi/2, -np.inf, logc)
        ratio = -np.expm1(q * logc) / s2
    ratio = np.where(s2 < 1e-16, q / 2, ratio)
    out = beta(q / 2, .5) * ratio / (2 * np.pi)
    return float(out) if out.ndim == 0 else out


def moment_observable(gamma, k=1):
    """Unbiased observable for E[(1-e^2)**(k-1/2)], integer k >= 1."""
    if not isinstance(k, (int, np.integer)) or k < 1:
        raise ValueError("k must be a positive integer")
    g = np.asarray(gamma, float)
    if np.any((g < 0) | (g > np.pi) | ~np.isfinite(g)):
        raise ValueError("gamma must be finite and in [0, pi]")
    return (np.pi / beta(k, .5) * np.sin(g) ** (2*k - 1)
            * (1 - (2*k - 1) * np.cos(g) ** 2))


def observable_bounds(k=1):
    """Exact minimum and maximum of the kth moment observable."""
    if not isinstance(k, (int, np.integer)) or k < 1:
        raise ValueError("k must be a positive integer")
    a = np.pi / beta(k, .5)
    b = 0 if k == 1 else 2 * ((2*k - 2)/(2*k + 1)) ** (k + .5)
    return -a*b, a


def estimate_moments(gamma, max_k=5, failure_probability=.05):
    """Means, standard errors and simultaneous distribution-free intervals.

    The Hoeffding intervals apply to independent observations of the ideal
    angle law. They are intersected with [0,1] only for display, never by
    clipping the point estimates. They do not account for model bias.
    """
    g = np.asarray(gamma, float)
    if g.ndim != 1 or len(g) < 2:
        raise ValueError("gamma must contain at least two observations")
    if not isinstance(max_k, (int, np.integer)) or max_k < 1:
        raise ValueError("max_k must be a positive integer")
    if not (0 < failure_probability < 1):
        raise ValueError("failure_probability must be in (0,1)")
    values = np.column_stack([moment_observable(g, k)
                              for k in range(1, max_k + 1)])
    means = values.mean(axis=0)
    cov = np.cov(values, rowvar=False, ddof=1) / len(g)
    cov = np.atleast_2d(cov)
    rows = []
    for idx, k in enumerate(range(1, max_k+1)):
        lo, hi = observable_bounds(k)
        radius = (hi-lo)*np.sqrt(np.log(2*max_k/failure_probability)/(2*len(g)))
        left, right = max(0., float(means[idx]-radius)), min(1., float(means[idx]+radius))
        rows.append(dict(k=k, power=2*k-1, estimate=float(means[idx]),
                         standard_error=float(np.sqrt(cov[idx, idx])),
                         simultaneous_hoeffding_radius=float(radius),
                         interval=[left, right] if left <= right else None,
                         empty_model_interval=bool(left > right)))
    return dict(n=len(g), failure_probability=failure_probability,
                moments=rows, covariance=cov.tolist())


def sample_exact(eccentricity, rng):
    """Fast sampler from the proven uniform/arcsine mixture factorisation."""
    e = np.asarray(eccentricity, float)
    if np.any((e < 0) | (e >= 1) | ~np.isfinite(e)):
        raise ValueError("eccentricity must be in [0,1)")
    t = np.sqrt((1-e)*(1+e)) * np.abs(np.sin(rng.uniform(0, 2*np.pi, e.shape)))
    u = rng.uniform(-1, 1, e.shape)
    return np.arctan2(t, u*np.sqrt(np.maximum(1-t*t, 0)))


def sample_kepler(eccentricity, rng, phase="mean", return_state=False):
    """Independent physical simulator; does not call the analytic kernel.

    Solves Kepler's equation, rotates orbital position and velocity, then
    measures their angle on the sky. `phase='true'` intentionally violates
    the theorem and is provided only for sensitivity experiments.
    """
    e = np.asarray(eccentricity, float)
    if np.any((e < 0) | (e >= 1) | ~np.isfinite(e)):
        raise ValueError("eccentricity must be in [0,1)")
    if phase == "mean":
        m = rng.uniform(0, 2*np.pi, e.shape)
        anom = np.where(e < .8, m, np.pi)
        for _ in range(80):
            step = (anom-e*np.sin(anom)-m)/(1-e*np.cos(anom))
            anom -= step
            if np.max(np.abs(step)) < 2e-14:
                break
        residual = float(np.max(np.abs(anom-e*np.sin(anom)-m)))
        if residual > 1e-11:
            raise ArithmeticError(f"Kepler solver residual {residual}")
    elif phase == "true":
        f = rng.uniform(0, 2*np.pi, e.shape)
        anom = np.arctan2(np.sqrt((1-e)*(1+e))*np.sin(f), e+np.cos(f))
        residual = None
    else:
        raise ValueError("phase must be 'mean' or 'true'")
    # Orientation: uniform argument of pericentre and uniform cos(i). The
    # longitude of the ascending node is omitted because a common rotation of
    # the sky plane leaves the projected angle unchanged; with it, the frame
    # would be exactly Haar-distributed on SO(3).
    j = np.sqrt((1-e)*(1+e))
    rx, ry = np.cos(anom)-e, j*np.sin(anom)
    vx, vy = -np.sin(anom), j*np.cos(anom)
    omega = rng.uniform(0, 2*np.pi, e.shape)
    ci = rng.uniform(-1, 1, e.shape)
    cw, sw = np.cos(omega), np.sin(omega)
    r1, r2 = rx*cw-ry*sw, ci*(rx*sw+ry*cw)
    v1, v2 = vx*cw-vy*sw, ci*(vx*sw+vy*cw)
    dot = r1*v1+r2*v2
    cross = np.abs(r1*v2-r2*v1)
    gamma = np.arctan2(cross, dot)
    if return_state:
        return dict(gamma=gamma, cos_inclination=ci,
                    projected_separation=np.hypot(r1,r2),
                    kepler_residual=residual)
    return gamma


def odd_sine_moment(eccentricity, k):
    """Exact E[sin(Gamma)**(2k+1)|e], k >= 1."""
    j = np.sqrt(1-np.asarray(eccentricity, float)**2)
    return beta(k,.5)/np.pi * sum(j**(2*r+1)/(2*r+1) for r in range(k))


def _check_gamma(gamma):
    g = np.asarray(gamma, float)
    if np.any((g < 0) | (g > np.pi) | ~np.isfinite(g)):
        raise ValueError("gamma must be finite and in [0, pi]")
    return g


def fractional_moment_observable(gamma, s):
    """Unbiased observable for E[(1-e^2)**(s/2)] = E[j**s], real s > -1.

    W_s(gamma) = pi/B((s+1)/2,1/2) |sin gamma|**s (1 - s cos(gamma)**2).
    It is bounded when s >= 0 and reduces to moment_observable(gamma, k)
    at s = 2k-1. For -1 < s < 0 it is unbounded near gamma = 0, pi.
    """
    s = float(s)
    if not np.isfinite(s) or s <= -1:
        raise ValueError("s must be finite and greater than -1")
    g = _check_gamma(gamma)
    return (np.pi / beta((s + 1) / 2, .5) * np.abs(np.sin(g)) ** s
            * (1 - s * np.cos(g) ** 2))


def fractional_observable_bounds(s):
    """Exact range [min, max] of fractional_moment_observable for s >= 0."""
    s = float(s)
    if not np.isfinite(s) or s < 0:
        raise ValueError("bounds are finite only for s >= 0")
    a = np.pi / beta((s + 1) / 2, .5)
    if s <= 1:
        return (a * (1 - s) if s == 0 else 0.0), a
    return -2 * a * ((s - 1) / (s + 2)) ** ((s + 2) / 2), a


def eccentricity_square_observable(gamma):
    """Unbiased bounded observable for E[e^2]; equals 1 + x - x^2, x = cos 2gamma.

    Equivalently 1/2 + cos(2 gamma) - cos(4 gamma)/2. Range [-1, 5/4].
    """
    x = np.cos(2 * _check_gamma(gamma))
    return 1 + x - x * x


def eccentricity_fourth_observable(gamma):
    """Unbiased bounded observable for E[e^4]; 1/3 + cos 2g - cos(6g)/3.

    Equivalently 1/3 + 2x - 4x^3/3 with x = cos 2gamma; range
    [1/3 - 2*sqrt(2)/3, 1/3 + 2*sqrt(2)/3].
    """
    x = np.cos(2 * _check_gamma(gamma))
    return 1 / 3 + 2 * x - 4 * x ** 3 / 3


def powerlaw_mean_circularity(alpha):
    """E[sqrt(1-e^2)] for f(e)=(alpha+1)e**alpha: (alpha+1)/2 B((alpha+1)/2,3/2)."""
    return (alpha + 1) / 2 * beta((alpha + 1) / 2, 1.5)


def powerlaw_alpha_from_sin3(gamma, upper=1e4):
    """Moment estimator of alpha solving E[j](alpha) = mean of (pi/2) sin^3.

    Returns nan when the sample mean lies outside the attainable range (0,1);
    that outcome is itself evidence against the ideal power-law model.
    """
    m = float(np.mean(moment_observable(_check_gamma(gamma), 1)))
    if not (powerlaw_mean_circularity(upper) < m < 1):
        return float("nan")
    return brentq(lambda a: powerlaw_mean_circularity(a) - m, -1 + 1e-12,
                  upper, xtol=1e-12)


def powerlaw_alpha_from_e2(gamma):
    """Closed-form moment estimator: E[e^2]=(alpha+1)/(alpha+3) inverted.

    Returns nan when the estimated mean square eccentricity is outside (0,1).
    """
    m = float(np.mean(eccentricity_square_observable(gamma)))
    if not (0 < m < 1):
        return float("nan")
    return (3 * m - 1) / (1 - m)


def powerlaw_fisher_information(alpha):
    """Fisher information per angle for alpha in the closed power-law family.

    Computed by adaptive quadrature in c = |cos gamma|; a floating-point
    value, not an interval certificate.
    """
    if not np.isfinite(alpha) or alpha <= -1:
        raise ValueError("alpha must be finite and greater than -1")
    q = alpha + 1
    const = 0.5 * (digamma(q / 2) - digamma(q / 2 + .5))

    def integrand(c):
        cq = c ** q
        score = const - cq * np.log(c) / (1 - cq)
        dens = beta(q / 2, .5) / np.pi * (1 - cq) / ((1 - c * c) ** 1.5)
        return score * score * dens
    return sum(quad(integrand, a, b, limit=400, epsabs=1e-14, epsrel=1e-10)[0]
               for a, b in [(0, .5), (.5, .9), (.9, 1)])


def abel_inverse_left_cdf(angle_density, r, singular_points=(), h=1e-6):
    """Recover F([0,r)) from an angle density by the inverse formula.

    Implements F([0,r)) = (2/pi) int_0^r t[-Y'(t)] / sqrt(r^2-t^2) dt with
    Y(c) = pi (1-c^2) p(arccos c), using t = r sin(theta) and a centred
    finite difference for Y'. `singular_points` lists eccentricity atoms, at
    which -Y' has an integrable inverse-square-root singularity. This is a
    numerical demonstration of the theorem, not a regularised estimator.
    """
    if not (0 < r < 1):
        raise ValueError("r must be in (0,1)")

    def Y(c):
        return np.pi * (1 - c * c) * angle_density(np.arccos(c))

    def integrand(th):
        t = r * np.sin(th)
        lo, hi = max(t - h, 0.0), min(t + h, 1.0)
        return (2 / np.pi) * t * -(Y(hi) - Y(lo)) / (hi - lo)
    pts = sorted(float(np.arcsin(e / r)) for e in singular_points if 0 < e < r)
    return quad(integrand, 0, np.pi / 2, points=pts or None, limit=600,
                epsabs=1e-10)[0]


def sample_kepler_physical(eccentricity, semimajor_axis, rng, phase="mean"):
    """Physical Kepler sample with sizes, in units with G*M_total = 1.

    Returns the projected angle, projected separation and projected relative
    speed. Orientation as in sample_kepler (full Haar up to a sky rotation).
    Used only for the selection experiments; does not call the kernel.
    """
    e = np.asarray(eccentricity, float)
    a = np.asarray(semimajor_axis, float)
    e, a = np.broadcast_arrays(e, a)
    if np.any((e < 0) | (e >= 1) | ~np.isfinite(e)) or np.any(a <= 0):
        raise ValueError("need 0 <= e < 1 and a > 0")
    if phase == "mean":
        m = rng.uniform(0, 2*np.pi, e.shape)
        anom = np.where(e < .8, m, np.pi)
        for _ in range(80):
            step = (anom-e*np.sin(anom)-m)/(1-e*np.cos(anom))
            anom = anom - step
            if np.max(np.abs(step)) < 2e-14:
                break
    elif phase == "eccentric":
        anom = rng.uniform(0, 2*np.pi, e.shape)
    else:
        raise ValueError("phase must be 'mean' or 'eccentric'")
    j = np.sqrt((1-e)*(1+e))
    rx, ry = a*(np.cos(anom)-e), a*j*np.sin(anom)
    speed = np.sqrt(1/a)/(1-e*np.cos(anom))
    vx, vy = -speed*np.sin(anom), speed*j*np.cos(anom)
    omega = rng.uniform(0, 2*np.pi, e.shape)
    ci = rng.uniform(-1, 1, e.shape)
    cw, sw = np.cos(omega), np.sin(omega)
    r1, r2 = rx*cw-ry*sw, ci*(rx*sw+ry*cw)
    v1, v2 = vx*cw-vy*sw, ci*(vx*sw+vy*cw)
    gamma = np.arctan2(np.abs(r1*v2-r2*v1), r1*v1+r2*v2)
    return dict(gamma=gamma, projected_separation=np.hypot(r1, r2),
                projected_speed=np.hypot(v1, v2))


def npmle_eccentricity(gamma, grid=None, iterations=5000, tol=1e-10):
    """Approximate grid maximum-likelihood fit of the eccentricity law (EM).

    Mixture weights on a fixed eccentricity grid (default 201 points on
    [0, 0.995]) are fitted by EM with the exact kernel. Returns a dict with the
    grid, weights, log-likelihood per observation at the returned weights,
    iterations used, and the optimality residual max_k mean(L_k / mix) - 1,
    which is <= 0 exactly at the grid maximum-likelihood solution. Kernel
    values are capped at 1e6 (only the circular grid point at 90 degrees can
    reach this); the number of capped entries is reported. Consistency and
    rates are not proved in this package.
    """
    g = _check_gamma(gamma)
    if grid is None:
        grid = np.linspace(0, .995, 201)
    grid = np.asarray(grid, float)
    with np.errstate(over="ignore"):
        L = angle_pdf(g[:, None], grid[None, :])
    capped = int(np.sum(L > 1e6))
    L = np.minimum(L, 1e6)
    w = np.full(len(grid), 1/len(grid))
    old = -np.inf
    used = iterations
    for it in range(iterations):
        mix = L @ w
        ll = float(np.mean(np.log(mix)))
        if ll - old < tol:
            used = it
            break
        old = ll
        w = w * (L.T @ (1 / mix)) / len(g)
        w = w / w.sum()
    mix = L @ w
    ll = float(np.mean(np.log(mix)))
    residual = float(np.max(L.T @ (1 / mix)) / len(g) - 1)
    return dict(grid=grid, weights=w, loglik_per_obs=ll, iterations=used,
                optimality_residual=residual, capped_kernel_entries=capped)


def step_cdf(grid, weights, r):
    """CDF F(r) = sum of weights at grid points <= r of a discrete law."""
    grid = np.asarray(grid, float); weights = np.asarray(weights, float)
    return np.array([weights[grid <= x].sum() for x in np.atleast_1d(r)])
