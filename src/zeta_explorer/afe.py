from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.zeta import _out, loggamma

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
