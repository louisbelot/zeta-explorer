from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.zeros import illinois, illinois_many

_TAIL = 40.0


def phi(u: ArrayLike, n_max: int = 30) -> np.ndarray:
    w = np.asarray(u)
    n2 = (np.arange(1, n_max + 1) ** 2).reshape((-1,) + (1,) * w.ndim)
    q = np.exp(-math.pi * n2 * np.exp(4 * w))
    return np.sum((2 * math.pi**2 * n2 * n2 * np.exp(9 * w) - 3 * math.pi * n2 * np.exp(5 * w)) * q,
                  axis=0)


@dataclass(frozen=True)
class _Contour:
    theta: float
    h: float
    u: np.ndarray
    G: np.ndarray


def _contour(t: float, x_max: float) -> _Contour:
    delta = math.pi / 8 * min(1.0, 4.0 / max(x_max, 1e-9))
    theta = math.pi / 8 - delta
    h = delta / 8
    s = math.sin(4 * delta)
    U = 0.25 * math.log(_TAIL / (math.pi * s)) + 0.25
    u = np.arange(0.0, U, h)
    n_max = math.ceil(math.sqrt(_TAIL / (math.pi * s))) + 2
    w = u + 1j * theta
    G = np.exp(t * w * w) * phi(w, n_max)
    G[0] /= 2
    return _Contour(theta, h, u, G)


def _evaluate(z: ArrayLike, t: float, order: int) -> complex | np.ndarray:
    zz = np.array(z, dtype=complex, ndmin=1).ravel()
    sign = np.where(zz.real < 0, -1.0, 1.0)
    zz = zz * sign
    c = _contour(t, float(np.max(np.abs(zz.real)) + np.max(np.abs(zz.imag))))
    e = np.exp(1j * np.outer(zz, c.u))
    up = c.G[None, :] * e * (1j * (c.u + 1j * c.theta)) ** order
    down = np.conj(c.G)[None, :] / e * (1j * (-c.u + 1j * c.theta)) ** order
    out = 0.5 * c.h * np.exp(-c.theta * zz) * np.sum(up + down, axis=1) * sign**order
    out = out.reshape(np.shape(z))
    return complex(out) if out.ndim == 0 else out


def ht(z: ArrayLike, t: float) -> complex | np.ndarray:
    return _evaluate(z, t, 0)


def ht_prime(z: ArrayLike, t: float) -> complex | np.ndarray:
    return _evaluate(z, t, 1)


def ht_second(z: ArrayLike, t: float) -> complex | np.ndarray:
    return _evaluate(z, t, 2)


def real_zeros(t: float, lo: float, hi: float, step: float = 0.05) -> np.ndarray:
    x = np.arange(lo, hi + step / 2, step)
    f = np.asarray(ht(x, t)).real
    sgn = np.sign(f)
    i = np.flatnonzero(sgn[:-1] * sgn[1:] < 0)
    if len(i) == 0:
        return np.empty(0)
    scale = np.abs(f).max()

    def g(q):
        return np.asarray(ht(q, t)).real / scale

    return illinois_many(g, x[i], x[i + 1], f[i] / scale, f[i + 1] / scale, rtol=1e-13)


def newton(z0: complex, t: float, tol: float = 1e-12, max_iter: int = 50) -> complex:
    z = complex(z0)
    for _ in range(max_iter):
        step = ht(z, t) / ht_prime(z, t)
        z -= step
        if abs(step) < tol * max(1.0, abs(z)):
            return z
    raise RuntimeError(f"Newton did not converge from {z0}")



def critical_point(t: float, lo: float, hi: float) -> float:
    fa, fb = (np.asarray(ht_prime([lo, hi], t)).real)
    scale = max(abs(fa), abs(fb))

    def g(q):
        return np.asarray(ht_prime(q, t)).real / scale

    return float(illinois_many(g, [lo], [hi], [fa / scale], [fb / scale], rtol=1e-14)[0])


def collision(lo: float, hi: float, t_lo: float, t_hi: float) -> tuple[float, float]:
    ref = abs(ht(lo, t_hi))

    def m(t):
        x = critical_point(t, lo, hi)
        return ht(x, t).real / ref

    t = illinois(m, t_lo, t_hi, m(t_lo), m(t_hi), tol=1e-13)
    return t, critical_point(t, lo, hi)

