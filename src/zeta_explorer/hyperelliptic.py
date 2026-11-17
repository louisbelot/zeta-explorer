from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from zeta_explorer.funcfield import _newton
from zeta_explorer.gf import Poly, character_sums, field, monics, squarefree
from zeta_explorer.twistmoments import _mul, _series


def ensemble(q: int, g: int) -> list[Poly]:
    return [D for D in monics(q, 2 * g + 1) if squarefree(D, q)]


@dataclass(frozen=True)
class Ensemble:
    q: int
    g: int
    D: list[Poly]
    c: np.ndarray

    @property
    def central(self) -> np.ndarray:
        g, w = self.g, np.full(self.g + 1, 2.0)
        w[g] = 1.0
        return self.c[:, : g + 1] @ (w * float(self.q) ** (-np.arange(g + 1) / 2))

    def angles(self) -> np.ndarray:
        q, g = self.q, self.g
        j = np.arange(g - 1, -1, -1)
        full = np.hstack([self.c, self.c[:, j] * float(q) ** (g - j)])
        out = np.empty((len(full), g))
        for i, c in enumerate(full):
            u = c / float(q) ** (np.arange(2 * g + 1) / 2)
            out[i] = np.sort(np.abs(np.angle(np.roots(u[::-1]))))[1::2]
        return out


def lpolys(q: int, g: int) -> Ensemble:
    Ds = ensemble(q, g)
    S = []
    for m in range(1, g + 1):
        F = field(q, m)
        S.append(character_sums(F, Ds, np.bincount(F.orbits[1])))
    S = np.array(S).T
    return Ensemble(q, g, Ds, np.array([_newton(list(map(int, s))) for s in S]))


def irreducible_counts(q: int, n: int) -> list[int]:
    mu = [0, 1] + [0] * n
    for k in range(2, n + 1):
        mu[k] = -sum(mu[d] for d in range(1, k) if k % d == 0)
    return [0] + [sum(mu[d] * q ** (k // d) for d in range(1, k + 1) if k % d == 0) // k
                  for k in range(1, n + 1)]


def recipe_first(q: int, g: int) -> float:
    n = 2 * g + 1
    G = _series(n)
    G[0, :] = [1, q] + [q**k - q ** (k - 1) for k in range(2, n + 1)]
    pi = irreducible_counts(q, g)
    for d in range(1, g // 2 + 1):
        local = _series(n)
        local[0, 0] = 1
        for k in range(1, g // (2 * d) + 1):
            for j in range(0, n + 1, d):
                local[2 * k * d, j] = (-1) ** (j // d)
        for _ in range(pi[d]):
            G = _mul(G, local)
    c = G[:, n]
    w = [2] * g + [1]
    return sum(w[j] * int(c[j]) * float(q) ** (-j / 2) for j in range(g + 1)) / int(G[0, n])


def andrade_keating(q: int, g: int, terms: int = 40) -> float:
    pi = irreducible_counts(q, terms)
    log_p = sum(n * math.log1p(-1 / ((q**d + 1) * q**d)) for d, n in enumerate(pi) if d)
    tail = sum(n * d / ((q**d + 1) * q**d - 1) for d, n in enumerate(pi) if d)
    return math.exp(log_p) / 2 * (2 * g + 2 + 4 * tail)


def usp_moment(N: int, k: float) -> float:
    lg = math.lgamma
    return math.exp(2 * N * k * math.log(2) + sum(
        lg(1 + N + j) + lg(0.5 + k + j) - lg(0.5 + j) - lg(1 + k + N + j) for j in range(1, N + 1)))
