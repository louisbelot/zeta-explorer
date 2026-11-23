from __future__ import annotations

import math
from dataclasses import dataclass

import mpmath
import numpy as np
from flint import arb, ctx, dirichlet_char

from zeta_explorer.selfdual import turing_bound

EPS = (0.1, 0.2, 0.3, 0.4, 0.5)


def characters(q: int) -> list[dirichlet_char]:
    out = []
    for n in range(1, q):
        if math.gcd(n, q) != 1 or pow(n, -1, q) < n:
            continue
        chi = dirichlet_char(q, n)
        if chi.is_primitive():
            out.append(chi)
    return out


def z_sign(chi: dirichlet_char, t: float, max_bits: int = 512) -> int:
    old = ctx.prec
    bits = 64
    try:
        while bits <= max_bits:
            ctx.prec = bits
            z = chi.hardy_z(arb(t)).real
            if z > 0:
                return 1
            if z < 0:
                return -1
            bits *= 2
        return 0
    finally:
        ctx.prec = old


def isolate(chi: dirichlet_char, lo: float, hi: float, h: float) -> list[tuple[float, float]]:
    n = max(2, math.ceil((hi - lo) / h) + 1)
    out = []
    prev_t, prev_s = None, 0
    for t in np.linspace(lo, hi, n):
        s = z_sign(chi, float(t))
        if s == 0:
            t = float(t) + (hi - lo) / (7 * n)
            s = z_sign(chi, t)
        if s and prev_s and s != prev_s:
            out.append((prev_t, float(t)))
        if s:
            prev_t, prev_s = float(t), s
    return out


def refine(chi: dirichlet_char, a: float, b: float, tol: float = 1e-12) -> tuple[float, float]:
    sa = z_sign(chi, a)
    while b - a > tol:
        m = (a + b) / 2
        sm = z_sign(chi, m)
        if sm == 0:
            return m, m
        if sm == sa:
            a = m
        else:
            b = m
    return a, b


def main_term(T: float, q: int, a: int) -> float:
    g = mpmath.loggamma(mpmath.mpc(0.25 + a / 2, T / 2)).imag
    return float(T / math.pi * math.log(q / math.pi) + 2 / math.pi * g)


def bmor(T: float, q: int, a: int) -> tuple[float, float, float]:
    ell = math.log(q * (T + 2) / (2 * math.pi))
    smooth = T / math.pi * math.log(q * T / (2 * math.pi * math.e)) - (1 - 2 * a) / 4
    return smooth, 0.22737 * ell + 2 * math.log1p(ell) - 0.5, ell / math.log(2 + ell)


@dataclass(frozen=True)
class Certificate:
    q: int
    n: int
    a: int
    t1: float
    t2: float
    zeros: list[tuple[float, float]]
    area: float
    bound: float

    @property
    def missing(self) -> float:
        return (self.bound - self.area) / (self.t2 - self.t1)

    @property
    def complete(self) -> bool:
        return self.missing < 1

    def _ends(self) -> tuple[np.ndarray, np.ndarray]:
        lo = np.array([min(abs(a), abs(b)) for a, b in self.zeros])
        hi = np.array([max(abs(a), abs(b)) for a, b in self.zeros])
        return lo, hi

    def straddles(self, T: float) -> bool:
        lo, hi = self._ends()
        return bool(((lo <= T) & (T < hi)).any())

    def count(self, T: float) -> int:
        if not self.complete or T > self.t1:
            raise ValueError("count is certified only up to t1")
        if self.straddles(T):
            raise ValueError(f"T = {T} falls inside an isolating interval")
        return int((self._ends()[1] <= T).sum())

    def gammas(self) -> np.ndarray:
        return np.array([(a + b) / 2 for a, b in self.zeros])


def certify(
    chi: dirichlet_char, T: float, window: float = 10.0, per_gap: int = 8, max_per_gap: int = 64
) -> Certificate:
    q, a = int(chi.modulus()), int(chi.parity())
    t1, t2 = float(T), float(T) + window
    bound = 2 * min(turing_bound(q, (a,), t1, t2, e) for e in EPS)
    while True:
        h = 2 * math.pi / (per_gap * math.log(q * t2 / (2 * math.pi)))
        zeros = isolate(chi, -t2, 0.0, h) + isolate(chi, 0.0, t2, h)
        area = window_area(zeros, q, a, t1, t2)
        cert = Certificate(q, int(chi.number()), a, t1, t2, zeros, area, bound)
        if cert.complete or per_gap >= max_per_gap:
            return cert
        per_gap *= 2


def window_area(zeros: list[tuple[float, float]], q: int, a: int, t1: float, t2: float) -> float:
    x, w = np.polynomial.legendre.leggauss(64)
    ts = (t1 + t2) / 2 + (t2 - t1) / 2 * x
    smooth = (t2 - t1) / 2 * float(w @ np.array([main_term(t, q, a) for t in ts]))
    u = np.array([max(abs(lo), abs(hi)) for lo, hi in zeros])
    u = u[u <= t2]
    return float(np.sum(t2 - np.maximum(u, t1))) - smooth
