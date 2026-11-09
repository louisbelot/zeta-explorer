from __future__ import annotations

import math

import mpmath
import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.zeros import illinois_many

_LOG_16 = math.log(16)
_LOG_PI = math.log(math.pi)
_HALF_LOG_2PI = 0.5 * math.log(2 * math.pi)
_CHUNK = 1 << 21


def log_m0(s: ArrayLike) -> np.ndarray:
    s = np.asarray(s, dtype=complex)
    return (np.log(s) + np.log(s - 1) - _LOG_16 - s / 2 * _LOG_PI + _HALF_LOG_2PI
            + (s / 2 - 0.5) * np.log(s / 2) - s / 2)


def alpha(s: ArrayLike) -> np.ndarray:
    s = np.asarray(s, dtype=complex)
    return 1 / (2 * s) + 1 / (s - 1) + 0.5 * np.log(s / (2 * math.pi))


def log_mt(s: ArrayLike, t: float) -> np.ndarray:
    return t / 4 * alpha(s) ** 2 + log_m0(s)


def log_bt(x: ArrayLike, y: float, t: float) -> np.ndarray:
    return log_mt((1 + y - 1j * np.asarray(x, dtype=float)) / 2, t)


def n_terms(x: ArrayLike, t: float) -> np.ndarray:
    return np.floor(np.sqrt(np.asarray(x, dtype=float) / (4 * math.pi) + t / 16)).astype(np.int64)


def _check(x: np.ndarray, y: float, t: float) -> None:
    if not 0 <= t <= 0.5 or not 0 <= y <= 1 or np.any(x < 200):
        raise ValueError("the effective bounds need 0 <= t <= 1/2, 0 <= y <= 1 and x >= 200")


def _by_terms(x: np.ndarray, t: float, block) -> np.ndarray:
    out = np.empty(len(x), dtype=complex)
    N = n_terms(x, t)
    for n in np.unique(N):
        idx = np.flatnonzero(N == n)
        ln = np.log(np.arange(1, n + 1, dtype=float))
        for part in np.array_split(idx, len(idx) * n // _CHUNK + 1):
            out[part] = block(part, ln)
    return out


def _shape(out: np.ndarray, x: ArrayLike):
    out = out.reshape(np.shape(x))
    return out[()] if out.ndim == 0 else out


def ft(x: ArrayLike, y: float, t: float) -> complex | np.ndarray:
    xs = np.atleast_1d(np.asarray(x, dtype=float)).ravel()
    _check(xs, y, t)
    s = (1 + y - 1j * xs) / 2
    ss = s + t / 2 * alpha(s)
    kappa = t / 2 * (alpha((1 - y + 1j * xs) / 2) - alpha((1 + y + 1j * xs) / 2))
    gamma = np.exp(log_mt((1 - y + 1j * xs) / 2, t) - log_mt(s, t))

    def block(i, ln):
        b = np.exp(t / 4 * ln * ln)
        a = np.exp(-np.outer(ss[i], ln)) @ b
        c = np.exp(-np.outer(np.conj(ss[i]) + kappa[i], ln)) @ (b * np.exp(y * ln))
        return a + gamma[i] * c

    return _shape(_by_terms(xs, t, block), x)


def c0(p: ArrayLike) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    d = 2 * np.cos(math.pi * p)
    num = np.exp(1j * math.pi * (p * p / 2 + 3 / 8)) - 1j * math.sqrt(2) * np.cos(math.pi * p / 2)
    near = np.abs(d) < 1e-8
    return np.where(near, (1 - 1j) / 4, num / np.where(near, 1.0, d))


def _log_m0_reduced(s: np.ndarray) -> np.ndarray:
    sigma, tau = s.real, s.imag
    return (np.log(np.abs(s)) + np.log(np.abs(s - 1)) - _LOG_16 - sigma / 2 * _LOG_PI
            + _HALF_LOG_2PI + (sigma / 2 - 0.5) * np.log(np.abs(s) / 2)
            + tau / 2 * np.arctan(sigma / tau) - sigma / 2)


def _c_scale(xs: np.ndarray, y: float, t: float) -> tuple[np.ndarray, np.ndarray]:
    T = xs / 2 + math.pi * t / 8
    s = (1 + y - 1j * xs) / 2
    d = (_log_m0_reduced(1j * T) - _log_m0_reduced(s) - t / 4 * (alpha(s) ** 2).real
         + t * math.pi**2 / 64 - t * math.pi**2 / 32)
    return T, np.exp(d)


def c_amplitude(x: ArrayLike, y: float, t: float) -> np.ndarray:
    xs = np.asarray(x, dtype=float)
    _check(np.atleast_1d(xs), y, t)
    T, scale = _c_scale(xs, y, t)
    a = np.sqrt(T / (2 * math.pi))
    N = np.floor(a)
    sign = np.where(N % 2 == 0, 1.0, -1.0)
    return 2 * sign * scale * np.imag(c0(1 - 2 * (a - N)) * np.exp(-1j * np.arctan(T)))


def c_ratio(x: ArrayLike, y: float, t: float) -> complex | np.ndarray:
    xs = np.asarray(x, dtype=float)
    return (np.exp(-1j * math.pi * y / 8) * c_amplitude(xs, y, t)
            * np.exp(-1j * log_bt(xs, y, t).imag))[()]


def ratio(x: ArrayLike, y: float, t: float, corrected: bool = True) -> complex | np.ndarray:
    out = ft(x, y, t)
    return out - c_ratio(x, y, t) if corrected else out


def error_ab(x: ArrayLike, y: float, t: float) -> float | np.ndarray:
    xs = np.atleast_1d(np.asarray(x, dtype=float)).ravel()
    _check(xs, y, t)
    s = (1 + y - 1j * xs) / 2
    sigma = (s + t / 2 * alpha(s)).real
    kappa = np.abs(t / 2 * (alpha((1 - y + 1j * xs) / 2) - alpha((1 + y + 1j * xs) / 2)))
    g = np.exp((log_mt((1 - y + 1j * xs) / 2, t) - log_mt(s, t)).real)

    def block(i, ln):
        n = len(ln)
        xi = xs[i, None]
        w = np.exp(t / 4 * ln * ln - np.outer(sigma[i], ln))
        lead = 1 + g[i, None] * n ** kappa[i, None] * np.exp(y * ln)
        tail = np.expm1((t * t / 16 * (np.log(xi / (4 * math.pi)) - 2 * ln) ** 2 + 0.626)
                        / (xi - 6.66))
        return np.sum(lead * w * tail, axis=1)

    return _shape(_by_terms(xs, t, block).real, x)


def error_c0(x: ArrayLike, y: float, t: float) -> float | np.ndarray:
    xs = np.asarray(x, dtype=float)
    _check(np.atleast_1d(xs), y, t)
    L = np.log(xs / (4 * math.pi))
    N = n_terms(xs, t)
    return np.exp(-(1 + y) / 4 * L - t / 16 * L * L + 1.24 * (3**y + 3**-y) / (N - 0.125)
                  + (3 * np.hypot(L, math.pi / 2) + 10.44) / (xs - 12))


def error_c(x: ArrayLike, y: float, t: float) -> float | np.ndarray:
    xs = np.asarray(x, dtype=float)
    _check(np.atleast_1d(xs), y, t)
    T, scale = _c_scale(xs, y, t)
    a = np.sqrt(T / (2 * math.pi))
    return scale * (1.24 * (3**y + 3**-y) / (a - 0.125) + 2 * 1.73 / (xs / 2 - 6))


def error(x: ArrayLike, y: float, t: float, corrected: bool = True) -> float | np.ndarray:
    tail = error_c(x, y, t) if corrected else error_c0(x, y, t)
    return error_ab(x, y, t) + tail


def _log_mt_mp(s, t):
    a = 1 / (2 * s) + 1 / (s - 1) + mpmath.log(s / (2 * mpmath.pi)) / 2
    m0 = (mpmath.log(s) + mpmath.log(s - 1) - mpmath.log(16) - s / 2 * mpmath.log(mpmath.pi)
          + mpmath.log(2 * mpmath.pi) / 2 + (s / 2 - mpmath.mpf(1) / 2) * mpmath.log(s / 2) - s / 2)
    return t / 4 * a * a + m0


def kernel_ratio(x: float, y: float, t: float, nodes: int = 48) -> complex:
    u, w = np.polynomial.hermite.hermgauss(nodes) if t > 0 else (np.zeros(1), np.ones(1))
    with mpmath.workdps(25 + int(math.log10(max(x, 10.0)))):
        X, Y, tm = mpmath.mpf(x), mpmath.mpf(y), mpmath.mpf(t)
        lb = _log_mt_mp((1 + Y - 1j * X) / 2, tm)
        acc = mpmath.mpc(0)
        for uk, wk in zip(u, w, strict=True):
            s = (1 - Y - 2 * mpmath.sqrt(tm) * mpmath.mpf(uk) + 1j * X) / 2
            lx = (mpmath.log(s) + mpmath.log(s - 1) - mpmath.log(16) - s / 2 * mpmath.log(mpmath.pi)
                  + mpmath.loggamma(s / 2) - lb)
            acc += mpmath.mpf(wk) * mpmath.exp(lx) * mpmath.zeta(s)
        return complex(acc / mpmath.sqrt(mpmath.pi) if t > 0 else acc)


def theta_t(x: ArrayLike, t: float) -> np.ndarray:
    return -log_bt(x, 0.0, t).imag


def theta_t_prime(x: ArrayLike, t: float) -> np.ndarray:
    s = (1 - 1j * np.asarray(x, dtype=float)) / 2
    a = alpha(s)
    da = -1 / (2 * s * s) - 1 / (s - 1) ** 2 + 1 / (2 * s)
    return (a + t / 2 * a * da).real / 2


def zt(x: ArrayLike, t: float, corrected: bool = True) -> float | np.ndarray:
    return np.real(np.exp(-1j * theta_t(x, t)) * ratio(x, 0.0, t, corrected))


def ideal_points(n: ArrayLike, t: float, tol: float = 1e-13) -> np.ndarray:
    target = (np.asarray(n, dtype=float) - 0.5) * math.pi
    x = 4 * math.pi * np.asarray(n, dtype=float) + 200
    for _ in range(100):
        dx = (theta_t(x, t) - target) / theta_t_prime(x, t)
        x = x - dx
        if np.all(np.abs(dx) < tol * x):
            return x
    raise RuntimeError("Newton didn't converge for some ideal points")


def displacement(x: ArrayLike, t: float) -> tuple[np.ndarray, np.ndarray]:
    k = theta_t(x, t) / math.pi + 0.5
    n = np.rint(k).astype(np.int64)
    return n, k - n


def real_zeros(t: float, lo: float, hi: float, per_gap: int = 8) -> np.ndarray:
    step = math.pi / float(theta_t_prime(hi, t)) / per_gap
    x = np.linspace(lo, hi, math.ceil((hi - lo) / step) + 1)
    f = zt(x, t)
    i = np.flatnonzero(f[:-1] * f[1:] < 0)
    if len(i) == 0:
        return np.empty(0)
    return illinois_many(lambda q: zt(q, t), x[i], x[i + 1], f[i], f[i + 1], rtol=1e-15)
