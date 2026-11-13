from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.hecke import Eigenform, eigenforms
from zeta_explorer.mellin import MellinLine
from zeta_explorer.selfdual import SelfDual, multiplicative
from zeta_explorer.sympower import _local, gamma_shifts


def sym_square(f: Eigenform) -> SelfDual:
    angle = dict(zip(f.p.tolist(), f.angles().tolist(), strict=True))

    def coeffs(M: int) -> np.ndarray:
        if M > f.P:
            raise ValueError(f"need lambda(p) up to {M}, have them up to {f.P}")
        return multiplicative(M, lambda p, K: _local(2, angle[p], K))

    return SelfDual(coeffs, 1.0, gamma_shifts(2, f.k), 1)


def harmonic_weight(f: Eigenform, L: SelfDual | None = None) -> float:
    L = sym_square(f) if L is None else L
    return 2 * math.pi**2 / ((f.k - 1) * L(1.0).real)


def unfold(L: SelfDual, t: ArrayLike) -> np.ndarray:
    return np.asarray(L.theta(np.asarray(t, dtype=float))) / math.pi


def low_zeros(L: SelfDual, x_max: float) -> np.ndarray:
    T = 1.0
    while unfold(L, T) < x_max + 0.5:
        T *= 1.5
    x = unfold(L, MellinLine(L, T).zeros(per_gap=32))
    return x[x <= x_max]


def count_zeros(L: SelfDual, T: float, sigma: float = 3.0, steps: int = 32) -> int:
    while True:
        v = L(np.linspace(sigma, 0.5, steps) + 1j * T)
        arg = np.unwrap(np.angle(v))
        if np.max(np.abs(np.diff(arg))) < 0.5:
            break
        if steps >= 1024:
            raise ArithmeticError(f"arg L(s) does not settle at t = {T}, too close to a zero")
        steps *= 2
    n = L.theta(T) / math.pi + arg[-1] / math.pi
    if abs(n - round(n)) > 0.1:
        raise ArithmeticError(f"N({T}) = {n:.3f} is not near an integer")
    return round(n)


@dataclass(frozen=True)
class Family:
    k: np.ndarray
    weight: np.ndarray
    zeros: list[np.ndarray]
    x_max: float

    @property
    def size(self) -> int:
        return len(self.k)

    def select(self, mask: ArrayLike) -> Family:
        m = np.asarray(mask, dtype=bool)
        zs = [z for z, keep in zip(self.zeros, m, strict=True) if keep]
        return Family(self.k[m], self.weight[m], zs, self.x_max)

    def shares(self, harmonic: bool = False) -> np.ndarray:
        w = self.weight if harmonic else np.ones(self.size)
        return w / w.sum()

    def _binned(self, bins: int, x_max: float,
                harmonic: bool) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        w = self.shares(harmonic) * self.size
        x = np.concatenate(self.zeros)
        wx = np.repeat(w, [len(z) for z in self.zeros])
        counts, edges = np.histogram(x, bins=bins, range=(0, x_max), weights=wx)
        inflation = np.sum(wx * wx) / np.sum(wx)
        return edges, counts, inflation

    def density(self, bins: int = 30, x_max: float = 3.0,
                harmonic: bool = False) -> tuple[np.ndarray, np.ndarray]:
        edges, counts, _ = self._binned(bins, x_max, harmonic)
        return (edges[:-1] + edges[1:]) / 2, counts / (self.size * (edges[1] - edges[0]))

    def chi2(self, model: Callable[[np.ndarray], np.ndarray], bins: int = 30, x_max: float = 3.0,
             harmonic: bool = False) -> float:
        edges, counts, inflation = self._binned(bins, x_max, harmonic)
        u, w = np.polynomial.legendre.leggauss(8)
        half = (edges[1] - edges[0]) / 2
        mid = (edges[:-1] + edges[1:]) / 2
        expect = self.size * half * (model(mid[:, None] + half * u) @ w)
        return float(np.sum((counts - expect) ** 2 / (inflation * expect)))

    def first(self) -> np.ndarray:
        return np.array([z[0] for z in self.zeros])

    def save(self, path: str) -> None:
        np.savez(path, k=self.k, weight=self.weight, zeros=np.concatenate(self.zeros),
                 lens=[len(z) for z in self.zeros], x_max=self.x_max)

    @classmethod
    def load(cls, path: str) -> Family:
        d = np.load(path)
        zs = np.split(d["zeros"], np.cumsum(d["lens"])[:-1])
        return cls(d["k"], d["weight"], zs, float(d["x_max"]))


def family(weights: ArrayLike, x_max: float = 3.0) -> Family:
    ks, ws, zs = [], [], []
    for k in np.asarray(weights, dtype=int).tolist():
        for f in eigenforms(k):
            L = sym_square(f)
            ks.append(k)
            ws.append(harmonic_weight(f, L))
            zs.append(low_zeros(L, x_max))
    return Family(np.array(ks), np.array(ws), zs, x_max)


def plancherel_mean(p: int, n: int) -> float:
    return p ** (-n / 2) if n % 2 == 0 else 0.0


def lambda_power(n: int, x: ArrayLike) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    prev, cur = np.zeros_like(x), np.ones_like(x)
    for _ in range(n):
        prev, cur = cur, x * cur - prev
    return cur

