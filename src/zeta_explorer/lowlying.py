from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import Character, is_fundamental, primitive_characters, quadratic
from zeta_explorer.lfunction import l_theta, l_zeros

_PSI = {
    0: -np.euler_gamma - math.pi / 2 - 3 * math.log(2),
    1: -np.euler_gamma + math.pi / 2 - 3 * math.log(2),
}


def unfold_l(gamma: ArrayLike, chi: Character) -> np.ndarray:
    th0 = l_theta(0.0, chi)
    return (np.asarray(l_theta(np.asarray(gamma, dtype=float), chi)) - th0) / math.pi


def low_zeros(chi: Character, x_max: float) -> np.ndarray:
    prim = chi.primitive()
    if math.log(prim.q / math.pi) + _PSI[prim.parity] <= 0:
        raise ValueError("theta is not increasing from t = 0, the conductor is too small")
    t = 1.0
    while unfold_l(t, prim) < x_max:
        t *= 2
    x = unfold_l(l_zeros(prim, 0.0, t), prim)
    return x[(x > 0) & (x <= x_max)]


def quadratic_family(lo: int, hi: int) -> list[Character]:
    ds = [d for m in range(lo, hi + 1) for d in (-m, m) if is_fundamental(d)]
    return [quadratic(d) for d in ds]


def prime_family(p: int) -> list[Character]:
    return [chi for _, chi in primitive_characters(p)]


@dataclass(frozen=True)
class OneLevel:
    zeros: list[np.ndarray]

    @property
    def x(self) -> np.ndarray:
        return np.concatenate(self.zeros)

    def density(self, bins: int = 30, x_max: float = 3.0) -> tuple[np.ndarray, np.ndarray]:
        counts, edges = np.histogram(self.x, bins=bins, range=(0, x_max))
        width = edges[1] - edges[0]
        return (edges[:-1] + edges[1:]) / 2, counts / (len(self.zeros) * width)

    def first(self) -> np.ndarray:
        return np.array([z[0] for z in self.zeros if len(z)])


def one_level(family: list[Character], x_max: float = 3.0) -> OneLevel:
    return OneLevel([low_zeros(chi, x_max) for chi in family])


def w_unitary(x: ArrayLike) -> np.ndarray:
    return np.ones_like(np.asarray(x, dtype=float))


def w_symplectic(x: ArrayLike) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    return 1 - np.sinc(2 * x)
