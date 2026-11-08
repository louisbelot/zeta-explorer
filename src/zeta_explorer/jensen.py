from __future__ import annotations

import math
from collections.abc import Iterable

import mpmath
from mpmath.calculus.quadrature import GaussLegendre


def _theta0(t: mpmath.mpf, eps: mpmath.mpf) -> mpmath.mpf:
    s = mpmath.mpf(0)
    k = 1
    while True:
        term = mpmath.exp(-mpmath.pi * k * k * t)
        s += term
        if term < eps * s:
            return s
        k += 1


def _upper(m_max: int, dps: int) -> float:
    target = math.log(10) * (dps + 20)
    U = 1.0
    while math.pi * math.exp(U) - m_max * math.log(U) - U / 4 < target:
        U += 0.25
    return U + 0.5


def theta_moments(ms: Iterable[int], dps: int = 60, panel: float = 0.05, degree: int = 4) -> list:
    ms = list(ms)
    with mpmath.workdps(dps + 10):
        eps = mpmath.mpf(10) ** (-dps - 10)
        rule = GaussLegendre(mpmath.mp).calc_nodes(degree, mpmath.mp.prec)
        U = _upper(max(ms), dps)
        panels = math.ceil(U / panel)
        h = mpmath.mpf(U) / panels
        u, g = [], []
        for p in range(panels):
            for x, w in rule:
                ui = h * (p + (x + 1) / 2)
                u.append(ui)
                g.append(w * h / 2 * mpmath.exp(ui / 4) * _theta0(mpmath.exp(ui), eps))
        return [mpmath.fsum(gi * ui**m for gi, ui in zip(g, u, strict=True)) for m in ms]


def taylor_gamma(ns: Iterable[int], dps: int = 60, **kw) -> list:
    ns = list(ns)
    ms = sorted({0} | {2 * n for n in ns} | {2 * n - 2 for n in ns if n > 0})
    F = dict(zip(ms, theta_moments(ms, dps, **kw), strict=True))
    out = []
    with mpmath.workdps(dps):
        for n in ns:
            if n == 0:
                out.append(4 - 2 * F[0])
                continue
            c = mpmath.factorial(n) / (mpmath.factorial(2 * n) * mpmath.mpf(2) ** (2 * n - 1))
            out.append(c * (32 * mpmath.binomial(2 * n, 2) * F[2 * n - 2] - F[2 * n]))
    return out
