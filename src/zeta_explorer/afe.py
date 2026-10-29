from __future__ import annotations

import cmath
import math

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import Character
from zeta_explorer.zeta import _out, loggamma

_LOSS = 3.0
_CUT = 40.0
_EPS = 2.0**-52
_TINY = 1e-300
_MAX_TERMS = 10_000


def gammaincc(z: complex, w: ArrayLike) -> complex | np.ndarray:
    ws = np.array(w, dtype=complex, ndmin=1).ravel()
    if z.real <= 0 or np.any((ws.real <= 0) & (ws.imag == 0)):
        raise ValueError("need Re z > 0 and w off the closed negative real axis")
    series, g = _kernel(z, ws)
    out = series + np.exp(z * np.log(ws) - ws - loggamma(z)) * g
    return _out(out.reshape(np.shape(w)), w)


def _kernel(z: complex, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    series = np.abs(w) < abs(z) + 1
    out = np.empty_like(w)
    out[series] = -_series(z, w[series]) / z
    out[~series] = _legendre(z, w[~series])
    return series, out


def _series(z: complex, w: np.ndarray) -> np.ndarray:
    term = np.ones_like(w)
    total = term.copy()
    for k in range(1, _MAX_TERMS):
        term = term * w / (z + k)
        total += term
        if np.all(np.abs(term) <= _EPS * np.abs(total)):
            return total
    raise RuntimeError("incomplete gamma series did not converge")


def _legendre(z: complex, w: np.ndarray) -> np.ndarray:
    b = w + 1 - z
    c = np.full_like(w, 1 / _TINY)
    d = 1 / b
    h = d.copy()
    out = np.empty_like(w)
    live = np.arange(len(w))
    for k in range(1, _MAX_TERMS):
        if len(live) == 0:
            return out
        a = -k * (k - z)
        b = b + 2
        d = 1 / (a * d + b)
        c = b + a / c
        delta = c * d
        h *= delta
        done = np.abs(delta - 1) <= 2 * _EPS
        out[live[done]] = h[done]
        keep = ~done
        live, b, c, d, h = live[keep], b[keep], c[keep], d[keep], h[keep]
    raise RuntimeError("Legendre's continued fraction did not converge")


def _angle(t: float) -> float:
    return math.pi / 2 if abs(t) * math.pi / 4 <= _LOSS else 2 * _LOSS / abs(t)


def afe_length(q: int, t: float) -> int:
    return math.isqrt(math.ceil((_CUT + _LOSS) * q / (math.pi * math.sin(_angle(t))))) + 1


def l_afe(s: ArrayLike, chi: Character) -> complex | np.ndarray:
    if chi.is_principal or not chi.is_primitive:
        raise ValueError("the smoothed sums need a primitive non-principal character")
    z = np.array(s, dtype=complex, ndmin=1).ravel()
    if np.any((z.real <= 0) | (z.real >= 1)):
        raise ValueError("l_afe works inside the critical strip 0 < Re s < 1")
    eps = chi.root_number()
    out = np.array([_afe_point(complex(si), chi, eps) for si in z], dtype=complex)
    return _out(out.reshape(np.shape(s)), s)


def _afe_point(s: complex, chi: Character, eps: complex) -> complex:
    q, a = chi.q, chi.parity
    phi = math.copysign(math.pi / 2 - _angle(s.imag), s.imag)
    n = np.arange(1, afe_length(q, s.imag) + 1)
    n = n[chi.num[n % q] >= 0]
    c = chi.values[n % q]
    w = math.pi / q * n.astype(float) ** 2 * cmath.exp(1j * phi)

    d1, m1 = _parts(s, a, n, c, w)
    d2, m2 = _parts(1 - s, a, n, c.conj(), w.conj())
    z1, z2 = (s + a) / 2, (1 - s + a) / 2
    lq = math.log(q / math.pi)
    lg1, lg2 = loggamma(z1), loggamma(z2)
    first = d1 + m1 * cmath.exp(-z1 * (lq - 1j * phi) - lg1)
    second = d2 * cmath.exp(lg2 - lg1) + m2 * cmath.exp(-z2 * (lq + 1j * phi) - lg1)
    return first + eps * cmath.exp((0.5 - s) * lq) * second


def _parts(
    s: complex, a: int, n: np.ndarray, c: np.ndarray, w: np.ndarray
) -> tuple[complex, complex]:
    series, g = _kernel((s + a) / 2, w)
    direct = np.sum(c[series] * np.exp(-s * np.log(n[series])))
    return complex(direct), complex(np.sum(c * n**a * np.exp(-w) * g))
