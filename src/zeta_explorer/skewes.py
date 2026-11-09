from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath
import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.multieval import zeros_fast
from zeta_explorer.riemann_siegel import z_rs
from zeta_explorer.spacing import first_zeros
from zeta_explorer.zeros import gram_points

BLOCK = 1 << 22


def lehman_sum(omega: ArrayLike, gamma: ArrayLike, alpha: float | None = None) -> np.ndarray:
    w = np.atleast_1d(np.asarray(omega, dtype=float))
    g = np.asarray(gamma, dtype=float)
    weight = np.ones_like(g) if alpha is None else np.exp(-(g**2) / (2 * alpha))
    c1 = weight / (0.25 + g**2)
    c2 = 2 * g * c1
    out = np.zeros_like(w)
    step = max(1, BLOCK // max(len(w), 1))
    for i in range(0, len(g), step):
        phase = np.outer(w, g[i : i + step])
        out -= np.cos(phase) @ c1[i : i + step] + np.sin(phase) @ c2[i : i + step]
    return out


def lehman_errors(
    omega: float, eta: float, alpha: float, T: float, A: float, improved: bool = True
) -> dict[str, mpmath.mpf]:
    lehman = 4 * A / omega <= alpha <= A**2 and 2 * A / alpha <= eta <= omega / 2
    if not (omega - eta > 1 and lehman):
        raise ValueError("parameters outside Lehman's conditions")
    if not 2 * math.pi * math.e < T <= A:
        raise ValueError("need 2 pi e < T <= A")
    if improved and omega - eta <= 437.5:
        raise ValueError("the improved s1 needs omega - eta > 437.5")
    with mpmath.workdps(30):
        lo, hi = mpmath.mpf(omega) - eta, mpmath.mpf(omega) + eta
        alpha, eta, T, A = (mpmath.mpf(v) for v in (alpha, eta, T, A))
        gauss = mpmath.exp(-alpha * eta**2 / 2)
        return {
            "s1": (mpmath.mpf("2.1457") if improved else mpmath.mpf("3.05")) / lo,
            "s2": 4 * hi * mpmath.exp(-lo / 6),
            "s3": 2 * gauss / (mpmath.sqrt(2 * mpmath.pi * alpha) * eta),
            "s4": mpmath.mpf("0.08") * mpmath.sqrt(alpha) * gauss,
            "s5": mpmath.exp(-(T**2) / (2 * alpha))
            * (alpha / (mpmath.pi * T**2) * mpmath.log(T / (2 * mpmath.pi)) + 8 * mpmath.log(T) / T
               + 4 * alpha / T**3),
            "s6": A * mpmath.log(A) * mpmath.exp(-(A**2) / (2 * alpha) + hi / 2)
            * (4 / mpmath.sqrt(alpha) + 15 * eta),
        }


def bulk_zeros(n: int, chunk: int = 20_000) -> np.ndarray:
    edges = np.arange(100, n + chunk + 100, chunk)
    gs = gram_points(edges)
    low = first_zeros(150)
    parts = [low[low < gs[0]]]
    for lo, hi, a, b in zip(edges[:-1], edges[1:], gs[:-1], gs[1:], strict=True):
        win, u = zeros_fast(int(lo) - 50, int(hi) + 200)
        u = u + win.t0
        parts.append(u[(u >= a) & (u < b)])
    return np.concatenate(parts)[:n]


def polish(gamma: ArrayLike, h: float = 1e-7, block: int = 1 << 16) -> np.ndarray:
    g = np.asarray(gamma, dtype=float)
    out = np.empty_like(g)
    for i in range(0, len(g), block):
        t = g[i : i + block]
        a, b = z_rs(t - h), z_rs(t + h)
        if np.any(np.signbit(a) == np.signbit(b)):
            raise ValueError(f"no sign change within {h} of some zero in block {i}")
        out[i : i + block] = t - h + 2 * h * a / (a - b)
    return out


def zero_error(omega: float, alpha: float, gamma: ArrayLike, delta: float) -> float:
    g = np.asarray(gamma, dtype=float) - delta
    q = 0.25 + g**2
    w = np.exp(-(g**2) / (2 * alpha))
    slope = ((g / alpha) * (1 + 2 * g) + omega + 2 + 2 * g * omega) / q + 2 * g * (1 + 2 * g) / q**2
    return float(delta * np.sum(w * slope))


def outward(x: mpmath.mpf, direction: int, digits: int = 9) -> str:
    e = int(mpmath.floor(mpmath.log10(x)))
    m = x / mpmath.mpf(10) ** (e - digits + 1)
    m = int(mpmath.floor(m) if direction < 0 else mpmath.ceil(m))
    s = str(m)
    return f"{s[0]}.{s[1:]}e+{e}"


@dataclass(frozen=True)
class Crossing:
    omega: float
    eta: float
    alpha: float
    T: float
    H: float
    dH: float
    errors: dict[str, mpmath.mpf]

    @property
    def R(self) -> float:
        return math.nextafter(float(sum(self.errors.values())), math.inf)

    @property
    def lower(self) -> float:
        return -1 + self.H - self.dH - self.R

    @property
    def interval(self) -> tuple[str, str]:
        with mpmath.workdps(30):
            w = mpmath.mpf(self.omega)
            return outward(mpmath.exp(w - self.eta), -1), outward(mpmath.exp(w + self.eta), 1)

    @property
    def excess(self) -> float:
        return self.lower * math.exp((self.omega - self.eta) / 2) / (self.omega + self.eta)


def crossing(
    omega: float,
    eta: float,
    alpha: float,
    gamma: ArrayLike,
    A: float,
    delta: float,
    improved: bool = True,
) -> Crossing:
    g = np.asarray(gamma, dtype=float)
    T = float(g[-1])
    H = float(lehman_sum(omega, g, alpha)[0])
    return Crossing(
        omega, eta, alpha, T, H, zero_error(omega, alpha, g, delta),
        lehman_errors(omega, eta, alpha, T, A, improved),
    )
