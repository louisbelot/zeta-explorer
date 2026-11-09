from __future__ import annotations

import functools

import numpy as np

from zeta_explorer.characters import factor
from zeta_explorer.elliptic import Curve, EllipticL
from zeta_explorer.explicit import primes
from zeta_explorer.selfdual import SelfDual, multiplicative, smoothed

__all__ = ["MU", "SymSquare", "coefficients", "smoothed"]

MU = (1.0, 1.0, 2.0)


@functools.cache
def coefficients(E: Curve, M: int) -> np.ndarray:
    N = E.conductor
    a = EllipticL(E).an(M).astype(float)

    def local(p: int, K: int) -> list[float]:
        x = a[p] ** 2 / p
        if N % p == 0:
            return [x**k for k in range(K + 1)]
        out = [1.0, x - 1]
        for _ in range(K - 1):
            prev2 = out[-3] if len(out) >= 3 else 0.0
            out.append((x - 1) * (out[-1] - out[-2]) + prev2)
        return out

    return multiplicative(M, local)


class SymSquare(SelfDual):
    def __init__(self, E: Curve):
        N = E.conductor
        if any(e > 1 for _, e in factor(N)):
            raise ValueError("sym^2 here needs a semistable curve")
        super().__init__(lambda M: coefficients(E, M), N * N, MU)
        self.E = E
        self.conductor = N * N

    def euler(self, s: complex, P: int = 1 << 16) -> complex:
        p = primes(P)
        x = EllipticL(self.E).an(P)[p].astype(float) ** 2 / p
        X = p.astype(float) ** -complex(s)
        bad = self.E.conductor % p == 0
        e = x - 1
        local = np.where(bad, 1 - x * X, 1 - e * X + e * X * X - X**3)
        return complex(np.exp(-np.sum(np.log(local))))
