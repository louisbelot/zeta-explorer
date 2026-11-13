from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from fractions import Fraction

import mpmath
import numpy as np
from flint import arb, arb_mat, fmpz_mat, fmpz_poly

from zeta_explorer.certified import precision
from zeta_explorer.explicit import primes


def dimension(k: int) -> int:
    if k < 12 or k % 2:
        return 0
    return k // 12 - (k % 12 == 2)


def _sigma(M: int, r: int) -> list[int]:
    s = [0] * (M + 1)
    for d in range(1, M + 1):
        dr = d**r
        for m in range(d, M + 1, d):
            s[m] += dr
    return s


@functools.cache
def _series(M: int) -> tuple[fmpz_poly, fmpz_poly, fmpz_poly]:
    E4 = fmpz_poly([1] + [240 * x for x in _sigma(M, 3)[1:]])
    E6 = fmpz_poly([1] + [-504 * x for x in _sigma(M, 5)[1:]])
    eta3 = [0] * (M + 1)
    n = 0
    while n * (n + 1) // 2 <= M:
        eta3[n * (n + 1) // 2] = (-1) ** n * (2 * n + 1)
        n += 1
    return E4, E6, fmpz_poly(eta3).pow_trunc(8, M).left_shift(1)


def basis(k: int, M: int) -> list[fmpz_poly]:
    E4, E6, D = _series(M)
    n = M + 1
    out = []
    for j in range(1, dimension(k) + 1):
        w = k - 12 * j
        b = w % 4 // 2
        a = (w - 6 * b) // 4
        out.append(D.pow_trunc(j, n).mul_low(E4.pow_trunc(a, n), n).mul_low(E6.pow_trunc(b, n), n))
    return out


def hecke_matrix(k: int, p: int = 2) -> fmpz_mat:
    d = dimension(k)
    g = basis(k, p * d)
    A = fmpz_mat(d, d, [g[j][i] for i in range(1, d + 1) for j in range(d)])
    B = fmpz_mat(d, d, [g[j][p * i] + (p ** (k - 1) * g[j][i // p] if i % p == 0 else 0)
                        for i in range(1, d + 1) for j in range(d)])
    return fmpz_mat([[int(x) for x in row] for row in A.solve(B).tolist()])


@dataclass(frozen=True)
class Eigenform:
    k: int
    a2: float
    p: np.ndarray
    lam: np.ndarray
    rad: float
    P: int

    def angles(self) -> np.ndarray:
        return np.arccos(np.clip(self.lam / 2, -1, 1))


def _eigenvectors(T: fmpz_mat, roots: list[arb]) -> list[arb_mat]:
    d = T.nrows()
    out = []
    for r in roots:
        S = arb_mat(T) - arb_mat(d, d, [r if i == j else 0 for i in range(d) for j in range(d)])
        if d == 1:
            out.append(arb_mat(1, 1, [1]))
            continue
        sub = arb_mat(d - 1, d - 1, [S[i, j] for i in range(1, d) for j in range(1, d)])
        rhs = arb_mat(d - 1, 1, [-S[i, 0] for i in range(1, d)])
        c = sub.solve(rhs)
        out.append(arb_mat(1, d, [arb(1)] + [c[i, 0] for i in range(d - 1)]))
    return out


@functools.cache
def eigenforms(k: int, P: int = 2048, tol: float = 1e-20) -> tuple[Eigenform, ...]:
    d = dimension(k)
    if d == 0:
        return ()
    ps = primes(P)
    g = basis(k, int(ps[-1]))
    T = hecke_matrix(k)
    poly = T.charpoly()
    G = fmpz_mat(d, len(ps), [g[i][p] for i in range(d) for p in ps.tolist()])
    bits = 128
    while True:
        with precision(bits):
            found = poly.complex_roots()
            if any(m > 1 for _, m in found):
                raise ValueError(f"T_2 has a repeated eigenvalue in weight {k}")
            roots = sorted((r.real for r, _ in found), key=lambda x: float(x.mid()))
            try:
                vectors = _eigenvectors(T, roots)
            except ZeroDivisionError:
                bits *= 2
                continue
            scale = [arb(p).sqrt() ** (k - 1) for p in ps.tolist()]
            coeffs = arb_mat(G)
            forms = []
            worst = 0.0
            for r, c in zip(roots, vectors, strict=True):
                a = c * coeffs
                lam = [a[0, j] / scale[j] for j in range(len(ps))]
                rad = max(float(x.rad()) for x in lam)
                worst = max(worst, rad)
                forms.append(Eigenform(k, float(r.mid()), ps,
                                       np.array([float(x.mid()) for x in lam]), rad, P))
        if worst < tol:
            return tuple(forms)
        bits *= 2


def kloosterman(m: int, n: int, c: int) -> float:
    return sum(math.cos(2 * math.pi * (m * x + n * pow(x, -1, c)) / c)
               for x in range(c) if math.gcd(x, c) == 1)


def petersson(k: int, m: int, n: int, C: int = 200, dps: int = 30) -> float:
    with mpmath.workdps(dps):
        x = 4 * mpmath.pi * mpmath.sqrt(m * n)
        s = mpmath.fsum(kloosterman(m, n, c) / c * mpmath.besselj(k - 1, x / c)
                        for c in range(1, C + 1))
        return float((m == n) + 2 * mpmath.pi * (-1) ** (k // 2) * s)


def hurwitz(N: int) -> Fraction:
    if N == 0:
        return Fraction(-1, 12)
    if N < 0 or N % 4 in (1, 2):
        return Fraction(0)
    h = Fraction(0)
    b = N % 2
    while 3 * b * b <= N:
        m = (b * b + N) // 4
        a = max(b, 1)
        while a * a <= m:
            if m % a == 0:
                c = m // a
                if a == b == c:
                    h += Fraction(1, 3)
                elif b == 0 and a == c:
                    h += Fraction(1, 2)
                elif b == 0 or a == b or a == c:
                    h += 1
                else:
                    h += 2
            a += 1
        b += 2
    return h


def _chebyshev(k: int, t: int, n: int) -> int:
    prev, cur = 0, 1
    for _ in range(k - 2):
        prev, cur = cur, t * cur - n * prev
    return cur


def trace(k: int, n: int) -> int:
    t_max = math.isqrt(4 * n)
    s = sum(_chebyshev(k, t, n) * hurwitz(4 * n - t * t) for t in range(-t_max, t_max + 1))
    tail = sum(min(d, n // d) ** (k - 1) for d in range(1, n + 1) if n % d == 0)
    total = -s / 2 - Fraction(tail, 2)
    if total.denominator != 1:
        raise ArithmeticError(f"trace of T_{n} on S_{k} came out as {total}")
    return int(total)
