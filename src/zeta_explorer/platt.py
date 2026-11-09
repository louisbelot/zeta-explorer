from __future__ import annotations

import math
import urllib.request
from dataclasses import dataclass

import numpy as np
from flint import arb

from zeta_explorer.certified import (
    Block,
    blocks_needed,
    gram_ball,
    gram_block,
    precision,
    z_ball,
    z_sign,
)
from zeta_explorer.riemann_siegel import theta
from zeta_explorer.zeros import gram_points

URL = "https://www.lmfdb.org/zeros/zeta/list?N={n}&limit={limit}"


def fetch(n: int, limit: int = 50) -> str:
    with urllib.request.urlopen(URL.format(n=n, limit=limit), timeout=60) as r:
        return r.read().decode()


def parse(text: str) -> tuple[list[int], list[str]]:
    rows = [line.split() for line in text.splitlines() if line.strip()]
    return [int(n) for n, _ in rows], [g for _, g in rows]


@dataclass(frozen=True)
class Isolation:
    gamma: arb
    root: arb
    sign: int

    @property
    def diff(self) -> float:
        return abs(float((self.root - self.gamma).mid()))


def isolate(g: str, delta: str = "1e-25", bits: int = 256) -> Isolation | None:
    with precision(bits):
        gamma, d = arb(g), arb(delta)
        lo, hi = z_ball(gamma - d, bits), z_ball(gamma + d, bits)
        if not ((lo > 0 and hi < 0) or (lo < 0 and hi > 0)):
            return None
        root = gamma - d + 2 * d * lo.mid() / (lo.mid() - hi.mid())
        return Isolation(gamma, root, 1 if lo > 0 else -1)


@dataclass(frozen=True)
class Window:
    a: int
    b: int
    before: int
    after: int
    blocks: list[Block]

    @property
    def changes(self) -> int:
        return sum(blk.changes for blk in self.blocks)

    @property
    def complete(self) -> bool:
        return self.changes == self.b - self.a

    def intervals(self) -> list[tuple[float, float]]:
        out = []
        for blk in self.blocks:
            pts = [(t, s) for t, s in zip(blk.points, blk.signs, strict=True) if s != 0]
            out += [(t0, t1) for (t0, s0), (t1, s1) in zip(pts[:-1], pts[1:], strict=True)
                    if s0 != s1]
        return out

    def index(self, t: float) -> int:
        for i, (t0, t1) in enumerate(self.intervals()):
            if t0 < t < t1:
                return self.a + 2 + i
        raise ValueError(f"no sign change around {t}")


def window(t_lo: float, t_hi: float, margin: int = 40, max_level: int = 8) -> Window:
    n0 = math.floor(theta(t_lo) / math.pi) - margin
    n1 = math.ceil(theta(t_hi) / math.pi) + margin
    ns = np.arange(n0, n1 + 1)
    g = [gram_ball(int(n), x) for n, x in zip(ns, gram_points(ns), strict=True)]
    gf = np.array([float(x.mid()) for x in g])
    sg = [z_sign(x) for x in g]
    if 0 in sg:
        raise RuntimeError("Z vanishes at a Gram point to working precision")
    good = [i for i in range(len(ns)) if sg[i] == (-1) ** int(ns[i] % 2)]
    blocks = [gram_block(j, k, gf, sg, max_level, n0)
              for j, k in zip(good[:-1], good[1:], strict=True)]
    rosser = [blk.changes >= blk.k - blk.j for blk in blocks]

    def run(lo: int, step: int) -> int:
        n = 0
        while 0 <= lo + step * n < len(blocks) and rosser[lo + step * n]:
            n += 1
        return n

    def below(i: int) -> bool:
        return run(i, -1) >= blocks_needed(gf[blocks[i].k - n0])

    def above(i: int) -> bool:
        r = run(i, 1)
        return r > 0 and r >= blocks_needed(gf[blocks[i + r - 1].k - n0])

    first = max(i for i, blk in enumerate(blocks) if gf[blk.k - n0] < t_lo)
    while not below(first):
        first -= 1
        if first < 0:
            raise RuntimeError("not enough Rosser blocks below the window")
    last = min(i for i, blk in enumerate(blocks) if gf[blk.j - n0] > t_hi)
    while not above(last):
        last += 1
        if last >= len(blocks):
            raise RuntimeError("not enough Rosser blocks above the window")
    inner = blocks[first + 1:last]
    return Window(blocks[first].k, blocks[last].j, run(first, -1), run(last, 1), inner)


@dataclass(frozen=True)
class PlattCheck:
    n: list[int]
    isolated: list[Isolation | None]
    window: Window

    @property
    def all_isolated(self) -> bool:
        return all(x is not None for x in self.isolated)

    @property
    def max_diff(self) -> float:
        return max(x.diff for x in self.isolated if x is not None)

    @property
    def indices_match(self) -> bool:
        return all(self.window.index(float(x.gamma.mid())) == n
                   for n, x in zip(self.n, self.isolated, strict=True) if x is not None)


def check(text: str, delta: str = "1e-25", bits: int = 256) -> PlattCheck:
    ns, gs = parse(text)
    iso = [isolate(g, delta, bits) for g in gs]
    w = window(float(gs[0]), float(gs[-1]))
    return PlattCheck(ns, iso, w)
