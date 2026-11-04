from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.explicit import primes
from zeta_explorer.zeta import digamma, loggamma, zeta

_CHUNK = 1 << 21


def zeta_logderiv(s: ArrayLike, r: float = 0.25, m: int = 32) -> complex | np.ndarray:
    z = np.array(s, dtype=complex, ndmin=1).ravel()
    if np.any(z == 1):
        raise ValueError("zeta'/zeta has a pole at s = 1")
    rot = np.exp(2j * math.pi * (np.arange(m) + 0.5) / m)
    ring = z[:, None] + r * rot
    g = (ring - 1) * zeta(ring.ravel()).reshape(ring.shape)
    dg = np.mean(g / rot, axis=1) / r
    out = dg / ((z - 1) * zeta(z)) - 1 / (z - 1)
    return complex(out[0]) if np.ndim(s) == 0 else out.reshape(np.shape(s))


def arithmetic(t: ArrayLike, p_max: int = 10**5) -> tuple[np.ndarray, np.ndarray]:
    ts = np.atleast_1d(np.asarray(t, dtype=float))
    p = primes(p_max).astype(float)
    lp = np.log(p)
    log_a = np.zeros(len(ts), dtype=complex)
    a_prime = np.zeros(len(ts), dtype=complex)
    rows = max(1, _CHUNK // len(p))
    for i in range(0, len(ts), rows):
        sl = slice(i, i + rows)
        u = np.exp(2j * ts[sl, None] * lp)
        log_a[sl] = np.log1p((1 - u) / (p * p - 1)).sum(axis=1)
        a_prime[sl] = (lp / ((p + 1) * (p * u - 1))).sum(axis=1)
    return np.exp(log_a), a_prime


def _theta(t: float, q: float, k: int) -> tuple[float, float]:
    w = (0.5 + k + 1j * t) / 2
    lq = math.log(q / math.pi)
    return t / 2 * lq + loggamma(w).imag, (lq + digamma(w).real) / 2


def _t_max(q: np.ndarray, kappa: np.ndarray, x_max: float) -> float:
    out = 0.0
    for k in np.unique(kappa):
        qk = float(q[kappa == k].min())
        t = math.pi * x_max / _theta(0.0, qk, k)[1]
        for _ in range(20):
            th, dth = _theta(t, qk, k)
            t -= (th - math.pi * x_max) / dth
        out = max(out, t)
    return 1.01 * out


def quadratic_density(
    ds: ArrayLike, x: ArrayLike, p_max: int = 10**4, n: int = 1500
) -> np.ndarray:
    ds = np.atleast_1d(np.asarray(ds, dtype=float))
    xs = np.asarray(x, dtype=float)
    q = np.abs(ds)
    kappa = (ds < 0).astype(int)
    if np.any(np.log(q / math.pi) + np.real(digamma((0.5 + kappa) / 2)) <= 0):
        raise ValueError("theta is not increasing from t = 0, a conductor is too small")
    t = np.linspace(0, _t_max(q, kappa, float(xs.max())), n + 1)[1:]
    a, a_prime = arithmetic(t, p_max)
    smooth = zeta_logderiv(1 + 2j * t) + a_prime
    twist = zeta(1 - 2j * t) * a
    out = np.zeros_like(xs)
    for k in (0, 1):
        w = (0.5 + k + 1j * t) / 2
        lg = np.imag(loggamma(w))
        dlg = np.real(digamma(w)) / 2
        for qi in q[kappa == k]:
            th = t / 2 * math.log(qi / math.pi) + lg
            b = smooth - np.exp(-2j * th) * twist
            rho = 1 + b.real / (math.log(qi / math.pi) / 2 + dlg)
            out += np.interp(xs, np.r_[0.0, th / math.pi], np.r_[0.0, rho])
    return out / len(ds)
