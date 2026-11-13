from __future__ import annotations

import functools
import math
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.explicit import primes
from zeta_explorer.hecke import Eigenform, dimension, eigenforms, petersson, trace
from zeta_explorer.mellin import MellinLine
from zeta_explorer.selfdual import SelfDual, _psi_gamma_factor, multiplicative
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


def jackson(x: ArrayLike, s: float) -> np.ndarray:
    return 0.75 * s * np.sinc(s * np.asarray(x, dtype=float) / 2) ** 4


def jackson_hat(xi: ArrayLike, s: float) -> np.ndarray:
    y = np.abs(2 * np.asarray(xi, dtype=float) / s)
    return np.where(y < 1, 1 - 1.5 * y**2 + 0.75 * y**3, np.where(y < 2, (2 - y) ** 3 / 4, 0.0))


Moments = Callable[[int, np.ndarray, int], np.ndarray]


def sato_tate(k: int, p: np.ndarray, j: int) -> np.ndarray:
    return np.full(len(p), float(j == 0))


def plancherel(k: int, p: np.ndarray, j: int) -> np.ndarray:
    return np.asarray(p, dtype=float) ** -j


def trace_means(k: int, p: np.ndarray, j: int) -> np.ndarray:
    return np.array([trace(k, q ** (2 * j)) / (q ** (j * (k - 1)) * dimension(k))
                     for q in np.asarray(p).tolist()])


def petersson_means(k: int, p: np.ndarray, j: int) -> np.ndarray:
    return np.array([petersson(k, 1, q ** (2 * j)) / petersson(k, 1, 1)
                     for q in np.asarray(p).tolist()])


@dataclass(frozen=True)
class Shape:
    k: int
    D: float
    t: np.ndarray
    x: np.ndarray
    rho: np.ndarray

    def g(self, t: ArrayLike, x0: ArrayLike, s: float) -> np.ndarray:
        u = self.D * np.asarray(t, dtype=float)[None, :]
        c = np.asarray(x0, dtype=float)[:, None]
        return self.D * (jackson(u - c, s) + jackson(u + c, s))


@functools.cache
def shape(k: int, x_end: float = 40.0, n: int = 40001) -> Shape:
    mu = gamma_shifts(2, k)
    L = SelfDual(lambda M: np.zeros(M + 1), 1.0, mu)
    T = 1.0
    while L.theta(T) / math.pi < x_end:
        T *= 1.5
    t = np.linspace(0, T, n)
    x = np.asarray(L.theta(t)) / math.pi
    rho = _psi_gamma_factor(0.5 + 1j * t, mu).real / math.pi
    D = 2 / float(np.interp(2.0, x, t))
    return Shape(k, D, t, x, rho)


def smoothed(F: Family, x0: ArrayLike, s: float, harmonic: bool = False) -> np.ndarray:
    x0 = np.atleast_1d(np.asarray(x0, dtype=float))
    share = F.shares(harmonic)
    out = np.zeros(len(x0))
    for k in np.unique(F.k).tolist():
        sh = shape(k)
        beyond = sh.x > F.x_max
        tail = np.trapezoid(sh.g(sh.t[beyond], x0, s) * sh.rho[beyond], sh.t[beyond], axis=1)
        for i in np.flatnonzero(F.k == k):
            t = np.interp(F.zeros[i], sh.x, sh.t)
            out += share[i] * (sh.g(t, x0, s).sum(axis=1) + tail)
    return out


def predicted(F: Family, x0: ArrayLike, s: float, means: Moments | None = None,
              harmonic: bool = False) -> np.ndarray:
    x0 = np.atleast_1d(np.asarray(x0, dtype=float))
    share = F.shares(harmonic)
    out = np.zeros(len(x0))
    for k in np.unique(F.k).tolist():
        sh = shape(k)
        if means is None:
            w = 1 - np.sinc(2 * sh.x)
            out += share[F.k == k].sum() * np.trapezoid(sh.g(sh.t, x0, s) * sh.rho * w, sh.t,
                                                        axis=1)
            continue
        value = np.trapezoid(sh.g(sh.t, x0, s) * sh.rho, sh.t, axis=1)
        U = 2 * math.pi * sh.D * s
        m = 1
        while m * math.log(2) < U:
            p = primes(int(math.exp(U / m)))
            lam = means(k, p, m) - means(k, p, m - 1) + 1
            u = m * np.log(p)
            ghat = 2 * np.cos(np.outer(x0, u) / sh.D) * jackson_hat(u / (2 * math.pi * sh.D), s)
            value -= ghat @ (np.log(p) * lam * p ** (-m / 2)) / (2 * math.pi)
            m += 1
        out += share[F.k == k].sum() * value
    return out
