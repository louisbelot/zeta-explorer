from __future__ import annotations

import cmath
import functools
import math

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.afe import gammaincc
from zeta_explorer.explicit import primes
from zeta_explorer.zeros import illinois_many
from zeta_explorer.zeta import loggamma

KAPPA = 5.5
A = 1 / (2 * math.pi)
_CUT = 56.0
_LOSS = 3.0


def _moduli(N: int) -> list[int]:
    small = primes(46341)
    out, M, p = [], 1, (1 << 31) - 1
    while M <= 4 * N**6:
        if all(p % q for q in small):
            out.append(p)
            M *= p
        p -= 2
    return out


def _residues(N: int, ps: list[int]) -> np.ndarray:
    k = np.arange(math.isqrt(2 * N) + 2)
    e = k * (k + 1) // 2
    c = np.where(k % 2, -1, 1) * (2 * k + 1)
    e, c = e[e < N], c[e < N]
    P = np.array(ps, dtype=np.int64)[:, None]
    pair = (e[:, None] + e[None, :]).ravel()
    keep = pair < N
    f = np.zeros(N, dtype=np.int64)
    np.add.at(f, pair[keep], np.outer(c, c).ravel()[keep])
    f = f[None, :] % P
    for _ in range(6):
        g = np.zeros_like(f)
        for s, ck in zip(e, c, strict=True):
            g[:, s:] += ck * f[:, : N - s]
        f = g % P
    return f


def _garner(r: np.ndarray, ps: list[int]) -> np.ndarray:
    v = []
    for j, p in enumerate(ps):
        acc = np.zeros(r.shape[1], dtype=np.int64)
        m = 1
        for i in range(j):
            acc = (acc + v[i] * m) % p
            m = m * ps[i] % p
        v.append((r[j] - acc) * pow(m, -1, p) % p)
    x = v[-1].astype(object)
    for i in range(len(ps) - 2, -1, -1):
        x = x * ps[i] + v[i].astype(object)
    M = math.prod(ps)
    return np.where(x > M // 2, x - M, x)


@functools.lru_cache(maxsize=4)
def _tau(N: int) -> np.ndarray:
    ps = _moduli(N)
    return np.concatenate([[0], _garner(_residues(N, ps), ps)]).astype(object)


def _size(N: int) -> int:
    return 1 << max(10, (N - 1).bit_length())


def tau(N: int) -> np.ndarray:
    if N < 1:
        raise ValueError("need N >= 1")
    return _tau(_size(N))[: N + 1].copy()


@functools.lru_cache(maxsize=4)
def _lam(N: int) -> np.ndarray:
    t = _tau(N)
    n = np.arange(1, N + 1, dtype=float)
    out = np.zeros(N + 1)
    out[1:] = np.array([float(x) for x in t[1:]]) * n**-KAPPA
    return out


def lam(N: int) -> np.ndarray:
    return _lam(_size(N))[: N + 1]


def theta(t: ArrayLike) -> float | np.ndarray:
    ts = np.asarray(t, dtype=float)
    out = ts * math.log(A) + np.imag(loggamma(KAPPA + 0.5 + 1j * ts))
    return float(out) if np.ndim(out) == 0 else out


def hardy_z(t: ArrayLike) -> float | np.ndarray:
    ts = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.array([_hardy_point(ti) for ti in ts])
    return float(out[0]) if np.ndim(t) == 0 else out.reshape(np.shape(t))


def _hardy_point(t: float) -> float:
    ang = math.pi / 2 if abs(t) * math.pi / 2 <= _LOSS else _LOSS / abs(t)
    y0 = cmath.exp(1j * math.copysign(math.pi / 2 - ang, t))
    M = math.ceil(_CUT * A / math.sin(ang)) + 1
    n = np.arange(1, M + 1)
    a = lam(M)[1:] / np.sqrt(n)
    S = np.sum(a * np.exp(-1j * t * np.log(n)) * gammaincc(complex(KAPPA + 0.5, t), n * y0 / A))
    return 2 * (cmath.exp(1j * theta(t)) * S).real


def completed(s: ArrayLike, y0: complex = 1.0) -> complex | np.ndarray:
    z = np.array(s, dtype=complex, ndmin=1).ravel()
    if np.any((z.real <= -KAPPA) | (z.real >= 1 + KAPPA)):
        raise ValueError("completed works for -11/2 < Re s < 13/2")
    y0 = complex(y0)
    if y0.real <= 0:
        raise ValueError("the split point needs Re y0 > 0")
    M = math.ceil(_CUT * A * max(abs(y0), 1 / abs(y0)) * abs(y0) / y0.real) + 1
    n = np.arange(1, M + 1)
    b = lam(M)[1:]
    log_r = np.log(A / n)
    out = np.empty(len(z), dtype=complex)
    for i, si in enumerate(z):
        a1, a2 = complex(si + KAPPA), complex(1 + KAPPA - si)
        g1 = cmath.exp(loggamma(a1)) * gammaincc(a1, n * y0 / A)
        g2 = cmath.exp(loggamma(a2)) * gammaincc(a2, n / (y0 * A))
        out[i] = np.sum(b * (np.exp(si * log_r) * g1 + np.exp((1 - si) * log_r) * g2))
    return complex(out[0]) if np.ndim(s) == 0 else out.reshape(np.shape(s))


def l_delta(s: ArrayLike) -> complex | np.ndarray:
    z = np.asarray(s, dtype=complex)
    return completed(z) / np.exp(z * math.log(A) + loggamma(z + KAPPA))


def zeros(T: float, per_gap: int = 8, rtol: float = 1e-14) -> np.ndarray:
    step = math.pi / max(1.0, math.log(T * A)) / per_gap
    ts = np.linspace(step, T, max(2, math.ceil(T / step)))
    z = hardy_z(ts)
    i = np.flatnonzero(z[:-1] * z[1:] < 0)
    return illinois_many(hardy_z, ts[i], ts[i + 1], z[i], z[i + 1], rtol)
