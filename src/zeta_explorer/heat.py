from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.zeros import illinois, illinois_many

_gauss = np.polynomial.legendre.leggauss

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


def _lost_pairs(before: np.ndarray, after: np.ndarray) -> list[tuple[float, float]]:
    kept = {int(np.argmin(np.abs(before - x))) for x in after}
    lost = [i for i in range(len(before)) if i not in kept]
    if len(lost) % 2:
        raise RuntimeError("a single real zero vanished, refine the t grid")
    return [(before[i], before[j]) for i, j in zip(lost[::2], lost[1::2], strict=True)]


def zero_paths(ts: ArrayLike, lo: float, hi: float, step: float = 0.01) -> list[np.ndarray]:
    ts = np.asarray(ts, dtype=float)
    if np.any(np.diff(ts) >= 0):
        raise ValueError("run the flow backwards: t must decrease")
    out, real, cplx = [], None, []
    for t in ts:
        r = real_zeros(t, lo, hi, step)
        cplx = [newton(z, t) for z in cplx]
        if real is not None and len(r) < len(real):
            for a, b in _lost_pairs(real, r):
                z = newton(complex((a + b) / 2, (b - a) / 2), t)
                if abs(z.imag) < 1e-8 * abs(z):
                    raise RuntimeError(f"no complex zero near {(a + b) / 2} at t = {t}")
                cplx.append(z.conjugate() if z.imag < 0 else z)
        out.append(np.concatenate([r.astype(complex), np.array(cplx, dtype=complex)]))
        real = r
    return out


def zero_velocity(x: ArrayLike, t: float) -> np.ndarray:
    return np.real(np.asarray(ht_second(x, t)) / np.asarray(ht_prime(x, t)))


def velocities(x: ArrayLike, which: ArrayLike, t: float = 0.0, nodes: int = 200) -> np.ndarray:
    x = np.sort(np.asarray(x, dtype=float))
    X = x[-1] + (x[-1] - x[-2]) / 2
    u, w = _gauss(nodes)
    s = (u + 1) / 2
    y = X / s
    rho = np.log(y / (4 * math.pi)) / (4 * math.pi) + t / (16 * y)
    out = []
    for k in np.atleast_1d(np.asarray(which, dtype=int)):
        xk = x[k]
        near = np.sum(1 / (xk - np.delete(x, k))) + np.sum(1 / (xk + x))
        tail = np.sum(w / 2 * rho * 2 * xk * X / (xk * xk * s * s - X * X))
        out.append(2 * (near + tail))
    return np.array(out)
