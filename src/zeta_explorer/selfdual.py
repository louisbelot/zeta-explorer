from __future__ import annotations

import functools
import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass

import mpmath
import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.explicit import primes
from zeta_explorer.zeros import illinois_many
from zeta_explorer.zeta import digamma, loggamma

_LOG_PI = math.log(math.pi)
_EPS = (0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5)


def _log_gamma_factor(s: np.ndarray, mu: tuple[float, ...]) -> np.ndarray:
    return sum(-(s + m) / 2 * _LOG_PI + loggamma((s + m) / 2) for m in mu)


def _psi_gamma_factor(s: np.ndarray, mu: tuple[float, ...]) -> np.ndarray:
    return sum(-_LOG_PI / 2 + digamma((s + m) / 2) / 2 for m in mu)


def smoothed(
    s: complex, a: np.ndarray, q: float, mu: tuple[float, ...], eps: int = 1,
    b: float = 30.0, h: float = 0.12, deriv: bool = False,
) -> complex | tuple[complex, complex]:
    c = max(s.real, 1 - s.real) + 0.5
    T = abs(s.imag)
    beta = math.copysign(sum(math.atan2(T, c - 0.5 + m) for m in mu) / 2, s.imag)
    k = math.ceil(math.sqrt(40 * b) / h)
    v = np.arange(-k, k + 1) * h
    u = c + 1j * v
    base = _log_gamma_factor(np.array([s]), mu)[0]
    log_g1 = u * u / b - 1j * beta * u + _log_gamma_factor(s + u, mu) - base
    log_g2 = u * u / b + 1j * beta * u + _log_gamma_factor(1 - s + u, mu) - base
    w1 = h / (2 * math.pi) * np.exp(log_g1 + u * math.log(q) / 2) / u
    w2 = h / (2 * math.pi) * np.exp(log_g2 + u * math.log(q) / 2) / u
    n = np.arange(1, len(a))
    ln = np.log(n)
    an = a[1:]
    e1 = np.exp(-np.outer(s + u, ln)) * an
    e2 = np.exp(-np.outer(1 - s + u, ln)) * an
    front = eps * q ** (0.5 - s)
    value = w1 @ e1.sum(axis=1) + front * (w2 @ e2.sum(axis=1))
    if not deriv:
        return complex(value)
    ps = _psi_gamma_factor(np.array([s]), mu)[0]
    d1 = w1 * (_psi_gamma_factor(s + u, mu) - ps)
    d2 = w2 * (-_psi_gamma_factor(1 - s + u, mu) - ps - math.log(q))
    dvalue = (d1 @ e1.sum(axis=1) - w1 @ (e1 @ ln)
              + front * (d2 @ e2.sum(axis=1) + w2 @ (e2 @ ln)))
    return complex(value), complex(dvalue)


def multiplicative(M: int, local: Callable[[int, int], Sequence[float]]) -> np.ndarray:
    out = np.ones(M + 1)
    out[0] = 0.0
    for p in primes(M).tolist():
        K = 1
        while p ** (K + 1) <= M:
            K += 1
        m = np.arange(p, M + 1, p)
        k = np.ones(len(m), dtype=np.int64)
        r = m // p
        while np.any(hit := r % p == 0):
            k[hit] += 1
            r[hit] //= p
        out[m] *= np.asarray(local(p, K))[k]
    return out


def _terms(s: complex, q: float, mu: tuple[float, ...]) -> int:
    size = q * math.prod(abs(s + m) / (2 * math.pi) + 1 for m in mu)
    return int(12 * math.sqrt(size)) + 64


class SelfDual:
    def __init__(
        self, coeffs: Callable[[int], np.ndarray], q: float, mu: tuple[float, ...],
        eps: int = 1, central: int | None = None,
    ):
        self.coeffs = coeffs
        self.q = float(q)
        self.mu = tuple(mu)
        self.eps = eps
        self.central = (1 - eps) // 2 if central is None else central

    def _a(self, s: complex) -> np.ndarray:
        M = _terms(s, self.q, self.mu)
        return self.coeffs(1 << max(10, (M - 1).bit_length()))[: M + 1]

    def __call__(self, s: ArrayLike) -> complex | np.ndarray:
        z = np.atleast_1d(np.asarray(s, dtype=complex))
        out = np.array([smoothed(si, self._a(si), self.q, self.mu, self.eps) for si in z.ravel()])
        return complex(out[0]) if np.ndim(s) == 0 else out.reshape(z.shape)

    def with_derivative(self, s: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        z = np.atleast_1d(np.asarray(s, dtype=complex)).ravel()
        pairs = [smoothed(si, self._a(si), self.q, self.mu, self.eps, deriv=True) for si in z]
        value, dvalue = np.array(pairs).T
        return value, dvalue

    def theta(self, t: ArrayLike) -> float | np.ndarray:
        ts = np.asarray(t, dtype=float)
        out = ts * math.log(self.q) / 2 + np.imag(_log_gamma_factor(0.5 + 1j * ts, self.mu))
        return float(out) if np.ndim(out) == 0 else out

    def phi(self, t: ArrayLike) -> float | np.ndarray:
        return self.theta(t) / math.pi - self.central / 2

    def density(self, t: float) -> float:
        psi = _psi_gamma_factor(np.array([0.5 + 1j * t]), self.mu)[0].real
        return (math.log(self.q) / 2 + psi) / math.pi

    def hardy(self, t: ArrayLike) -> float | np.ndarray:
        ts = np.atleast_1d(np.asarray(t, dtype=float))
        w = np.exp(1j * self.theta(ts)) * self(0.5 + 1j * ts)
        out = w.real if self.eps == 1 else w.imag
        return float(out[0]) if np.ndim(t) == 0 else out.reshape(np.shape(t))

    def zeros(self, T: float, per_gap: int = 8, rtol: float = 1e-12) -> np.ndarray:
        return scan_zeros(self.hardy, self.density, T, per_gap, rtol)


def scan_zeros(
    hardy: Callable[[np.ndarray], np.ndarray], density: Callable[[float], float], T: float,
    per_gap: int = 8, rtol: float = 1e-12,
) -> np.ndarray:
    step = 1 / (max(density(T), 0.25) * per_gap)
    ts = np.linspace(1e-3, T, max(2, math.ceil(T / step) + 1))
    z = hardy(ts)
    i = np.flatnonzero(z[:-1] * z[1:] < 0)
    return illinois_many(hardy, ts[i], ts[i + 1], z[i], z[i + 1], rtol)


def _a_eps(eps: float) -> float:
    return (0.5 + eps) * math.log(1 + 1 / (0.5 + eps)) + math.log(1.5 + eps)


@functools.cache
def c0(eps: float) -> float:
    zeta = mpmath.zeta

    def log_z(x):
        return mpmath.log(zeta(2 * x) / zeta(x))

    def log_big(x):
        return mpmath.log(zeta(x))

    a = 1 + mpmath.mpf(eps)
    slope = 2 * zeta(2 * a, derivative=1) / zeta(2 * a) - zeta(a, derivative=1) / zeta(a)
    out = (
        (a - 0.5) * log_big(a)
        + mpmath.quad(log_big, [a, 4, 12, mpmath.inf])
        - mpmath.quad(log_z, [a, a + 1])
        - mpmath.quad(log_z, [1.5, 4, 12, mpmath.inf])
        + _a_eps(eps) * slope
    )
    return float(out)


def s_integral_bound(L: SelfDual, t1: float, t2: float, eps: float = 0.5) -> float:
    if not 0 < eps <= 0.5:
        raise ValueError("eps must lie in (0, 1/2]")
    x2 = min(t1 * t1 - (2 + eps + m) ** 2 for m in L.mu)
    if x2 <= 25 or t2 < t1:
        raise ValueError("t1 is too small for these gamma shifts, or t2 < t1")
    r = len(L.mu)
    a = _a_eps(eps)

    def log_q(s: complex) -> float:
        return math.log(L.q) + sum(math.log(abs(s + m) / (2 * math.pi)) for m in L.mu)

    total = (
        (1 / 16 + eps * (1 + eps) / 4) * log_q(complex(3 + eps, t2))
        + (a - 1) / 2 * log_q(complex(1 + eps, t1))
        + r * c0(eps)
        + r * (8.3 + 0.09 * a + (0.5 + eps) * (2 + eps) / 5) / math.sqrt(x2)
    )
    return total / math.pi


@dataclass(frozen=True)
class Turing:
    t1: float
    t2: float
    found: int
    expected: float
    area: float
    bound: float

    @property
    def missing(self) -> float:
        return (self.bound - self.area) / (self.t2 - self.t1)

    @property
    def complete(self) -> bool:
        return self.missing < 1


def turing(L: SelfDual, zeros: ArrayLike, t1: float, t2: float) -> Turing:
    g = np.sort(np.asarray(zeros, dtype=float))
    g = g[g <= t2]
    x, w = np.polynomial.legendre.leggauss(64)
    t = (t1 + t2) / 2 + (t2 - t1) / 2 * x
    area = float(np.sum(t2 - np.maximum(g, t1)) - (t2 - t1) / 2 * (w @ L.phi(t)))
    bound = min(s_integral_bound(L, t1, t2, e) for e in _EPS)
    return Turing(t1, t2, int(np.sum(g <= t1)), float(L.phi(t1)), area, bound)

