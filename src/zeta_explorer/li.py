from __future__ import annotations

import math

import mpmath
import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.moments import EULER_GAMMA
from zeta_explorer.riemann_siegel import theta, theta_prime

_CHUNK = 1 << 22


def li_trend(n: ArrayLike) -> np.ndarray:
    n = np.asarray(n, dtype=float)
    return n / 2 * (np.log(n) - 1 + EULER_GAMMA - math.log(2 * math.pi))


def li_exact(n_max: int, dps: int = 40, nodes: int = 160, r: float = 0.5) -> np.ndarray:
    with mpmath.workdps(dps):
        logs = []
        for j in range(nodes):
            z = r * mpmath.expjpi(mpmath.mpf(2 * j) / nodes)
            s = 1 / (1 - z)
            xi = s * (s - 1) / 2 * mpmath.pi ** (-s / 2) * mpmath.gamma(s / 2) * mpmath.zeta(s)
            w = mpmath.log(xi)
            if logs:
                w += 2j * mpmath.pi * mpmath.nint((logs[-1] - w).imag / (2 * mpmath.pi))
            logs.append(w)
        out = []
        for n in range(1, n_max + 1):
            roots = (mpmath.expjpi(mpmath.mpf(-2 * j * n) / nodes) for j in range(nodes))
            c = mpmath.fdot(logs, roots)
            out.append(float((n * c / nodes / mpmath.mpf(r) ** n).real))
        return np.array(out)


def _pair_terms(n: np.ndarray, gamma: np.ndarray) -> np.ndarray:
    out = np.zeros(len(n))
    phi = np.arctan(0.5 / gamma)
    rows = max(1, _CHUNK // len(n))
    for i in range(0, len(phi), rows):
        out += (4 * np.sin(np.multiply.outer(phi[i : i + rows], n)) ** 2).sum(axis=0)
    return out


def li_tail(n: ArrayLike, T: float, nodes: int = 200) -> np.ndarray:
    n = np.atleast_1d(np.asarray(n, dtype=float))
    y, w = np.polynomial.legendre.leggauss(nodes)
    y, w = (y + 1) / 2, w / 2
    phi_T = math.atan(0.5 / T)
    phi = phi_T * y**4
    t = 0.5 / np.tan(phi)
    dens = theta_prime(t) / math.pi * 2 / np.sin(phi) ** 2 * 4 * phi_T * y**3 * w
    return (np.sin(np.multiply.outer(n, phi)) ** 2) @ dens


def li_coefficients(n: ArrayLike, gamma: ArrayLike) -> np.ndarray:
    n = np.atleast_1d(np.asarray(n, dtype=float))
    gamma = np.sort(np.asarray(gamma, dtype=float))
    T = float(gamma[-1])
    S = len(gamma) - float(theta(T)) / math.pi - 1
    edge = 4 * np.sin(n * math.atan(0.5 / T)) ** 2
    return _pair_terms(n, gamma) + li_tail(n, T) - edge * S
