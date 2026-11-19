from __future__ import annotations

from fractions import Fraction

from zeta_explorer.funcfield import LPoly, Twist


def power_sums(c: list[int] | tuple[int, ...], n: int) -> list[int]:
    c = list(c) + [0] * n
    p: list[int] = []
    for m in range(1, n + 1):
        p.append(-m * c[m] - sum(c[k] * p[m - k - 1] for k in range(1, m)))
    return p


def frob_traces(L: LPoly, n: int) -> list[Fraction]:
    return [Fraction(x, L.q**m) for m, x in enumerate(power_sums(L.c, n), 1)]


def family_traces(fam: list[Twist], n: int, eps: int | None = None) -> list[Fraction]:
    rows = [frob_traces(t.L, n) for t in fam if eps is None or t.L.eps == eps]
    return [sum(col, Fraction(0)) / len(rows) for col in zip(*rows, strict=True)]


def rmt_trace(dim: int, det: int, n: int) -> int:
    if n < 1:
        raise ValueError("n must be >= 1")
    even, low = n % 2 == 0, n <= dim - 2
    if dim % 2 == 0:
        if det == 1:
            return int(even and low)
        return 2 - low if even else 0
    v = 1 - (not even and low)
    return v if det == 1 else (-1) ** n * v
