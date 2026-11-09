from __future__ import annotations

import math
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass

import flint
import numpy as np
from flint import acb, acb_series, arb

from zeta_explorer.riemann_siegel import theta
from zeta_explorer.zeros import gram_points

TRUDGIAN_A = 2.067
TRUDGIAN_B = 0.059
T_MIN = 168 * math.pi
G_MINUS_1 = 9.666908056130192


@contextmanager
def precision(bits: int, cap: int | None = None) -> Iterator[None]:
    old = flint.ctx.prec, flint.ctx.cap
    flint.ctx.prec = bits
    if cap is not None:
        flint.ctx.cap = cap
    try:
        yield
    finally:
        flint.ctx.prec, flint.ctx.cap = old


def _arb(t: float | arb) -> arb:
    return t if isinstance(t, arb) else arb(t)


def theta_ball(t: float | arb, bits: int = 64) -> arb:
    with precision(bits):
        t = _arb(t)
        return acb(0.25, t / 2).lgamma().imag - t / 2 * arb.pi().log()


def theta_prime_ball(t: float | arb, bits: int = 64) -> arb:
    with precision(bits):
        t = _arb(t)
        return (acb(0.25, t / 2).digamma().real - arb.pi().log()) / 2


def z_ball(t: float | arb, bits: int = 64) -> arb:
    with precision(bits):
        t = _arb(t)
        return (acb(0, theta_ball(t, bits)).exp() * acb(0.5, t).zeta()).real


def zeta_ball(s: complex, bits: int = 64) -> acb:
    with precision(bits):
        return acb(s.real, s.imag).zeta()


def z_sign(t: float | arb, bits: int = 64, max_bits: int = 1024) -> int:
    while bits <= max_bits:
        z = z_ball(t, bits)
        if z > 0:
            return 1
        if z < 0:
            return -1
        bits *= 2
    return 0


def gram_ball(n: int, guess: float, bits: int = 128) -> arb:
    with precision(bits):
        target = n * arb.pi()
        m = arb(guess)
        for _ in range(8):
            m = (m - (theta_ball(m, bits) - target) / theta_prime_ball(m, bits)).mid()
        box = arb(m, arb(2) ** (16 - bits) * (1 + abs(m)))
        k = m - (theta_ball(m, bits) - target) / theta_prime_ball(box, bits)
        if not box.contains_interior(k):
            raise RuntimeError(f"could not isolate Gram point {n}")
        return k


def blocks_needed(g: float) -> int:
    L = math.log(g)
    return math.ceil(L * (TRUDGIAN_A + TRUDGIAN_B * L) / (6 * math.pi))


def count_sign_changes(signs: list[int]) -> int:
    s = [x for x in signs if x != 0]
    return sum(a != b for a, b in zip(s[:-1], s[1:], strict=True))


@dataclass(frozen=True)
class Block:
    j: int
    k: int
    points: list[float]
    signs: list[int]

    @property
    def changes(self) -> int:
        return count_sign_changes(self.signs)


@dataclass(frozen=True)
class Certificate:
    n: int
    g_n: arb
    blocks: list[Block]
    turing_blocks: int

    @property
    def zeros(self) -> int:
        return sum(b.changes for b in self.blocks)

    def log(self) -> dict:
        return {
            "n": self.n,
            "g_n": self.g_n.str(30, radius=True),
            "zeros_on_line": self.zeros,
            "turing_blocks": self.turing_blocks,
            "blocks": [
                {"j": b.j, "k": b.k, "points": b.points, "signs": b.signs} for b in self.blocks
            ],
        }


def gram_block(
    j: int, k: int, gf: np.ndarray, sg: list[int], max_level: int, first: int = -1
) -> Block:
    need = k - j
    pieces = max(need, 1)
    for _ in range(max_level):
        pieces *= 2
        inner = np.linspace(gf[j], gf[k], pieces + 1)[1:-1].tolist()
        signs = [sg[j], *(z_sign(t) for t in inner), sg[k]]
        if count_sign_changes(signs) >= need:
            break
    return Block(j + first, k + first, [float(gf[j]), *inner, float(gf[k])], signs)


def certify(T: float, max_level: int = 8) -> Certificate:
    T = max(T, T_MIN)
    n0 = math.ceil(theta(T) / math.pi)
    margin = 64
    while True:
        ns = range(-1, n0 + margin)
        guess = np.concatenate(([G_MINUS_1], gram_points(np.arange(n0 + margin))))
        g = [gram_ball(n, x) for n, x in zip(ns, guess, strict=True)]
        gf = np.array([float(x.mid()) for x in g])
        sg = [z_sign(x) for x in g]
        if 0 in sg:
            raise RuntimeError("Z vanishes at a Gram point to working precision")
        good = [i for i, n in enumerate(ns) if sg[i] == (-1) ** n]
        if good[0] != 0:
            raise RuntimeError("g_-1 is a bad Gram point")
        found = _turing(good, gf, sg, n0 + 1, max_level)
        if found is not None:
            break
        margin *= 2
    p, used = found
    pairs = [(j, k) for j, k in zip(good[:-1], good[1:], strict=True) if k <= p]
    cert = Certificate(p - 1, g[p], [gram_block(j, k, gf, sg, max_level) for j, k in pairs], used)
    if cert.zeros != p:
        raise RuntimeError(f"found {cert.zeros} sign changes below g_{p - 1}, expected {p}")
    return cert


def _turing(
    good: list[int], gf: np.ndarray, sg: list[int], i0: int, max_level: int
) -> tuple[int, int] | None:
    start, run = None, 0
    for j, k in zip(good[:-1], good[1:], strict=True):
        if j < i0:
            continue
        if start is None:
            start = j
        if k - j > 1 and gram_block(j, k, gf, sg, max_level).changes < k - j:
            start, run = None, 0
            continue
        run += 1
        if run >= blocks_needed(float(gf[k])):
            return start, run
    return None


def _dh_series(s: acb) -> tuple[acb, acb]:
    k = ((arb(10) - 2 * arb(5).sqrt()).sqrt() - 2) / (arb(5).sqrt() - 1)
    x = acb_series([s, 1])
    f = sum(
        (c * x.zeta(acb(a) / 5) for a, c in zip((1, 2, 3, 4), (1, k, -k, -1), strict=True)),
        acb_series([0]),
    )
    f = f * (-x * arb(5).log()).exp()
    c = f.coeffs()
    return c[0], c[1]


def dh_ball(s: complex, bits: int = 64) -> acb:
    with precision(bits, cap=2):
        return _dh_series(acb(s.real, s.imag))[0]


def krawczyk(s: complex, r: float, bits: int = 128) -> acb | None:
    with precision(bits, cap=2):
        m = acb(s.real, s.imag)
        box = acb(arb(s.real, r), arb(s.imag, r))
        fm, dm = _dh_series(m)
        y = 1 / dm.mid()
        _, d_box = _dh_series(box)
        k = m - y * fm + (1 - y * d_box) * (box - m)
        return k if box.contains_interior(k) else None
