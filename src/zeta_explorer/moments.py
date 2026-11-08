from __future__ import annotations

import math
from collections.abc import Callable
from fractions import Fraction
from functools import cache

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.explicit import ei, primes
from zeta_explorer.multieval import ZWindow
from zeta_explorer.riemann_siegel import TWO_PI
from zeta_explorer.zeta import bernoulli, zeta

EULER_GAMMA = 0.5772156649015329
LOW = 1000.0

CFKRS = {
    2: (
        0.0506605918211688857219397316048638,
        0.69886988487897996984709628427658502,
        2.425962198846682004756575310160663,
        3.227907964901254764380689851274668,
        1.312424385961669226168440066229978,
    ),
    3: (
        0.000005708527034652788398376841445252313,
        0.00040502133088411440331215332025984,
        0.011072455215246998350410400826667,
        0.14840073080150272680851401518774,
        1.0459251779054883439385323798059,
        3.984385094823534724747964073429,
        8.60731914578120675614834763629,
        10.274330830703446134183009522,
        6.59391302064975810465713392,
        0.9165155076378930590178543,
    ),
}


def moment_polynomial(k: int) -> np.ndarray:
    if k == 1:
        return np.array([1.0, 2 * EULER_GAMMA])
    if k in CFKRS:
        return np.array(CFKRS[k])
    raise ValueError("P_k is tabulated for k = 1, 2, 3")


def barnes_g(n: int) -> int:
    return math.prod(math.factorial(j) for j in range(n - 1))


def prime_square_tail(P: float) -> float:
    return float(-ei(-math.log(P)).real)


def arithmetic_factor(k: int, P: int = 10**6) -> float:
    x = 1.0 / primes(P)
    term = np.ones_like(x)
    total = np.ones_like(x)
    m = 0
    while np.max(term) > 1e-18 * np.min(total):
        term *= ((m + k) / (m + 1)) ** 2 * x
        total += term
        m += 1
    log_a = np.sum(k * k * np.log1p(-x) + np.log(total))
    tail = k * k * (k - 1) ** 2 / 4 * prime_square_tail(P)
    return math.exp(log_a - tail)


def leading_coefficient(k: int, P: int = 10**6) -> float:
    return arithmetic_factor(k, P) * barnes_g(k + 1) ** 2 / barnes_g(2 * k + 1)


def keating_snaith(k: int, x: ArrayLike, P: int = 10**6) -> np.ndarray:
    return leading_coefficient(k, P) * np.asarray(x, dtype=float) ** (k * k)


def _log_powers_integral(j: int, t: np.ndarray) -> np.ndarray:
    L = np.log(t / TWO_PI)
    return t * sum((-1) ** i * math.perm(j, i) * L ** (j - i) for i in range(j + 1))


def polynomial_integral(k: int, C: float, D: float) -> float:
    c = moment_polynomial(k)[::-1]
    out = 0.0
    for j, cj in enumerate(c):
        lo = 0.0 if C == 0 else float(_log_powers_integral(j, np.float64(C)))
        out += cj * (float(_log_powers_integral(j, np.float64(D))) - lo)
    return out


def _series_exp(c: np.ndarray) -> np.ndarray:
    out = np.zeros_like(c)
    out[0] = math.exp(c[0])
    for n in range(1, len(c)):
        out[n] = sum(j * c[j] * out[n - j] for j in range(1, n + 1)) / n
    return out


def smoothed_polynomial_integral(k: int, T: float) -> float:
    p = moment_polynomial(k)[::-1]
    deg = len(p) - 1
    c = np.zeros(deg + 1)
    c[1] = math.log(T / TWO_PI) - EULER_GAMMA
    for m in range(2, deg + 1):
        c[m] = (-1) ** m * float(zeta(m).real) / m
    e = _series_exp(c)
    return T * sum(p[j] * math.factorial(j) * e[j] for j in range(deg + 1))


@cache
def gregory_weights(m: int) -> tuple[float, ...]:
    A = [[Fraction(j) ** d for j in range(m)] for d in range(m)]
    b = [bernoulli(d + 1) / (d + 1) if d % 2 else Fraction(0) for d in range(m)]
    c = _solve(A, b)
    c[0] -= Fraction(1, 2)
    return tuple(float(1 + x) for x in c)


def _solve(A: list[list[Fraction]], b: list[Fraction]) -> list[Fraction]:
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(r for r in range(col, n) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        for r in range(n):
            if r != col and M[r][col] != 0:
                f = M[r][col] / M[col][col]
                M[r] = [a - f * v for a, v in zip(M[r], M[col], strict=True)]
    return [M[i][n] / M[i][i] for i in range(n)]


def z_power_integrals(
    C: float,
    D: float,
    ks: ArrayLike = (1, 2, 3),
    weight: Callable[[np.ndarray], np.ndarray] | None = None,
    h: float | None = None,
    window: float = 20_000.0,
    order: int = 8,
) -> np.ndarray:
    ks = np.atleast_1d(np.asarray(ks, dtype=int))
    if h is None:
        h = 1 / (2 * ks.max() * math.log(max(D, LOW) / TWO_PI))
    n = max(2 * order, math.ceil((D - C) / h))
    h = (D - C) / n
    g = np.array(gregory_weights(order))
    out = np.zeros(len(ks))
    edges = [C, LOW] if C < LOW < D else [C]
    edges += list(np.arange(max(C, LOW) + window, D, window)) + [D]
    for a, b in zip(edges[:-1], edges[1:], strict=True):
        i0 = math.ceil((a - C) / h - 1e-9)
        i1 = n if b == D else math.ceil((b - C) / h - 1e-9) - 1
        i = np.arange(i0, i1 + 1)
        t = C + h * i
        w = np.full(len(i), h)
        head, tail = i < order, i > n - order
        w[head] *= g[i[head]]
        w[tail] *= g[n - i[tail]]
        if weight is not None:
            w *= weight(t)
        if b <= LOW:
            z2 = np.abs(zeta(0.5 + 1j * t)) ** 2
        else:
            win = ZWindow(t[0] - 1, t[-1] + 1)
            z2 = win.z(t - win.t0) ** 2
        out += [np.dot(w, z2**k) for k in ks]
    return out


def smoothed_moments(T: float, ks: ArrayLike = (1, 2, 3), cut: float = 60.0, **kw) -> np.ndarray:
    return z_power_integrals(0.0, cut * T, ks, weight=lambda t: np.exp(-t / T), **kw)
