from __future__ import annotations

import functools
import math

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import kronecker
from zeta_explorer.elliptic import Curve, EllipticL, t_of_x
from zeta_explorer.explicit import primes
from zeta_explorer.symsquare import SymSquare
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


def _twist_locals(E: Curve, P: int) -> tuple[np.ndarray, np.ndarray, int, float]:
    p = primes(P)
    a = EllipticL(E).an(P)
    M = E.conductor
    good = p[p != M]
    return good.astype(float), a[good].astype(float) ** 2 / good, M, a[M] / math.sqrt(M)


def twist_arithmetic(
    E: Curve, t: ArrayLike, eta: int, P: int = 10**5
) -> tuple[np.ndarray, np.ndarray]:
    ts = np.atleast_1d(np.asarray(t, dtype=complex))
    p, lam2, M, lam_m = _twist_locals(E, P)
    lp = np.log(p)
    lam = np.sqrt(lam2)
    e = lam2 - 1
    a1 = np.zeros(len(ts), dtype=complex)
    log_a = np.zeros(len(ts), dtype=complex)
    rows = max(1, _CHUNK // len(p))
    for i in range(0, len(ts), rows):
        sl = slice(i, i + rows)
        phase = np.exp(1j * ts[sl, None] * lp)
        X = phase ** -2 / p
        Q = 1 - e * X + e * X * X - X**3
        g = p / (p + 1) * X * (2 + 2 * X - lam2) / ((1 + X) ** 2 - lam2 * X)
        a1[sl] = (lp * (g - X / (1 - X) + X * (e - 2 * e * X + 3 * X * X) / Q)).sum(axis=1)
        x = phase / np.sqrt(p)
        y = 1 / (phase * np.sqrt(p))
        even = 0.5 * ((1 - lam * y + y * y) / (1 - lam * x + x * x)
                      + (1 + lam * y + y * y) / (1 + lam * x + x * x))
        local = 1 + p / (p + 1) * (even - 1)
        Y = x * x
        fix = (1 - y * y) * (1 - e * Y + e * Y * Y - Y**3) / (
            (1 - 1 / p) * (1 - e / p + e / p**2 - 1 / p**3)
        )
        log_a[sl] = np.log(local * fix).sum(axis=1)
    beta = lam_m * lam_m
    lm = math.log(M)
    xm = np.exp(-(0.5 + 1j * ts) * lm)
    Xm = xm * xm
    a1 += lm * (-eta * lam_m * xm / (1 - eta * lam_m * xm) - Xm / (1 - Xm)
                + beta * Xm / (1 - beta * Xm))
    ym = np.exp(-(0.5 - 1j * ts) * lm)
    log_a += np.log((1 - eta * lam_m * xm) / (1 - eta * lam_m * ym)
                    * (1 - Xm) * (1 - beta * ym * ym)
                    / ((1 - 1 / M) * (1 - beta / M)))
    return a1, np.exp(log_a)


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


@functools.cache
def twist_r1(E: Curve, eta: int, P: int = 10**5) -> float:
    a1, _ = twist_arithmetic(E, 0.0, eta, P)
    value, dvalue = SymSquare(E).with_derivative(1.0)
    return float(np.euler_gamma - (dvalue[0] / value[0]).real - a1[0].real)


def twist_effective_n(L: EllipticL, P: int = 10**5) -> float:
    r1 = twist_r1(L.E, kronecker(L.d, L.E.conductor), P)
    return (math.log(L.A) - np.euler_gamma) / (2 * r1)


def twist_profile(
    E: Curve, t: ArrayLike, etas: tuple[int, ...] = (1, -1), P: int = 10**5
) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    t = np.asarray(t, dtype=float)
    s = 1 + 2j * t
    L2 = SymSquare(E)
    sym, dsym = L2.with_derivative(s)
    smooth = -zeta_logderiv(s) + dsym / sym
    twist = zeta(s) * np.conj(sym) / L2(1.0)
    out = {}
    for eta in etas:
        a1, a = twist_arithmetic(E, t, eta, P)
        out[eta] = smooth + a1, twist * a
    return out


def twist_rho(
    t: np.ndarray, b1: np.ndarray, k: np.ndarray, log_a: float, root_number: int
) -> tuple[np.ndarray, np.ndarray]:
    th = t * log_a + np.imag(loggamma(1 + 1j * t))
    b = b1 - root_number * np.exp(-2j * th) * k
    return th / math.pi, 1 + b.real / (log_a + np.real(digamma(1 + 1j * t)))


def twist_density(
    E: Curve, ds: ArrayLike, x: ArrayLike, P: int = 10**5, n: int = 800
) -> np.ndarray:
    Ls = [EllipticL(E, int(d)) for d in np.atleast_1d(ds)]
    xs = np.asarray(x, dtype=float)
    t_max = 1.01 * t_of_x(float(xs.max()), min(Ls, key=lambda L: L.conductor))
    t = np.linspace(0, t_max, n + 1)[1:]
    etas = {L: kronecker(L.d, E.conductor) for L in Ls}
    profile = twist_profile(E, t, tuple(set(etas.values())), P)
    out = np.zeros_like(xs)
    for L, eta in etas.items():
        x, rho = twist_rho(t, *profile[eta], math.log(L.A), L.root_number)
        out += np.interp(xs, x, rho)
    return out / len(Ls)
