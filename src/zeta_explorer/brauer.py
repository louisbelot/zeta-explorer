from __future__ import annotations

import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np

from zeta_explorer.bsd import torsion
from zeta_explorer.funcfield import Twist, family, legendre
from zeta_explorer.gf import Poly
from zeta_explorer.heights import lstar, model


@dataclass(frozen=True)
class Ratio:
    D: Poly
    q: int
    chi: int
    rank: int
    lstar: Fraction
    tamagawa: int
    sha_reg: Fraction

    @property
    def log_h(self) -> float:
        return self.chi * math.log(self.q)

    @property
    def bs(self) -> float:
        return math.log(self.sha_reg) / self.log_h


def ratio(t: Twist) -> Ratio:
    q = t.L.q
    M = model(t.D, q)
    tam = math.prod(f.c for f in M.fibers.values())
    ls = lstar(t.L)
    sr = ls * Fraction(q) ** (M.chi - 1) * torsion(t.D, q) ** 2 / tam
    return Ratio(t.D, q, M.chi, t.rank, ls, tam, sr)


@dataclass(frozen=True)
class Slice:
    q: int
    n: int
    chi: int
    bs: np.ndarray
    sha_reg: np.ndarray
    rank: np.ndarray
    log_lstar: np.ndarray
    log_tamagawa: np.ndarray

    @property
    def log_h(self) -> float:
        return self.chi * math.log(self.q)

    @property
    def trivial(self) -> float:
        return float(np.mean(self.sha_reg[self.rank == 0] == 1))


def scan(q: int, n: int) -> Slice:
    rs = [ratio(t) for t in family(legendre(q), n)]
    chi = {r.chi for r in rs}
    if len(chi) != 1:
        raise ArithmeticError("height varies inside one degree")
    return Slice(q, n, chi.pop(), np.array([r.bs for r in rs]),
                 np.array([float(r.sha_reg) for r in rs]), np.array([r.rank for r in rs]),
                 np.array([math.log(r.lstar) for r in rs]),
                 np.array([math.log(r.tamagawa) for r in rs]))
