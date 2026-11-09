from __future__ import annotations

import functools
import math

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import factor
from zeta_explorer.elliptic import Curve, EllipticL
from zeta_explorer.explicit import primes
from zeta_explorer.zeta import digamma, loggamma

MU = (1.0, 1.0, 2.0)
_LOG_PI = math.log(math.pi)


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
    beta = math.copysign(max(0.0, math.pi * len(mu) / 4 - 2 / T), s.imag) if T > 0 else 0.0
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


@functools.cache
def coefficients(E: Curve, M: int) -> np.ndarray:
    N = E.conductor
    a = EllipticL(E).an(M).astype(float)
    out = np.ones(M + 1)
    out[0] = 0.0
    for p in primes(M).tolist():
        x = a[p] ** 2 / p
        local = [1.0, x if N % p == 0 else x - 1]
        while p ** len(local) <= M:
            if N % p == 0:
                local.append(local[-1] * x)
            else:
                prev2 = local[-3] if len(local) >= 3 else 0.0
                local.append((x - 1) * (local[-1] - local[-2]) + prev2)
        m = np.arange(p, M + 1, p)
        k = np.ones(len(m), dtype=np.int64)
        r = m // p
        while np.any(hit := r % p == 0):
            k[hit] += 1
            r[hit] //= p
        out[m] *= np.array(local)[k]
    return out


def _terms(s: complex, q: float, mu: tuple[float, ...]) -> int:
    size = q * math.prod(abs(s + m) / (2 * math.pi) + 1 for m in mu)
    return int(12 * math.sqrt(size)) + 64


class SymSquare:
    def __init__(self, E: Curve):
        N = E.conductor
        if any(e > 1 for _, e in factor(N)):
            raise ValueError("sym^2 here needs a semistable curve")
        self.E = E
        self.conductor = N * N

    def _a(self, s: complex) -> np.ndarray:
        M = _terms(s, self.conductor, MU)
        return coefficients(self.E, 1 << max(10, (M - 1).bit_length()))[: M + 1]

    def __call__(self, s: ArrayLike) -> complex | np.ndarray:
        z = np.atleast_1d(np.asarray(s, dtype=complex))
        out = np.array([smoothed(si, self._a(si), self.conductor, MU) for si in z.ravel()])
        return complex(out[0]) if np.ndim(s) == 0 else out.reshape(z.shape)

    def with_derivative(self, s: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
        z = np.atleast_1d(np.asarray(s, dtype=complex)).ravel()
        pairs = [smoothed(si, self._a(si), self.conductor, MU, deriv=True) for si in z]
        value, dvalue = np.array(pairs).T
        return value, dvalue

    def euler(self, s: complex, P: int = 1 << 16) -> complex:
        p = primes(P)
        x = EllipticL(self.E).an(P)[p].astype(float) ** 2 / p
        X = p.astype(float) ** -complex(s)
        bad = self.E.conductor % p == 0
        e = x - 1
        local = np.where(bad, 1 - x * X, 1 - e * X + e * X * X - X**3)
        return complex(np.exp(-np.sum(np.log(local))))
