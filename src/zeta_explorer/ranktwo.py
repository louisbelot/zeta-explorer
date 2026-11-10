from __future__ import annotations

import math
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np

from zeta_explorer.elliptic import (
    Curve,
    EllipticL,
    central,
    hardy_z,
    low_zeros,
    second_central,
    twists,
)
from zeta_explorer.rmt import (
    Kernel,
    char_poly_at_one,
    conditioned_n,
    fit_n,
    gap_chi2,
    gap_ks,
    independent_n,
    interaction_n,
    jacobi_angles,
)

KERNELS = {"independent": independent_n, "interaction": interaction_n,
           "conditioned": conditioned_n}


@dataclass(frozen=True)
class Twist:
    d: int
    log_A: float
    second: float
    x: np.ndarray


def rank_two(L: EllipticL, x_max: float = 3.0) -> Twist:
    z = low_zeros(L, x_max)
    if z.rank != 2:
        raise ValueError(f"{L} has central order {z.rank}, not 2")
    return Twist(L.d, math.log(L.A), second_central(L), z.x)


def even_values(E: Curve, lo: int, hi: int) -> tuple[list[EllipticL], np.ndarray]:
    even = [L for L in twists(E, lo, hi, sign=1) if L.root_number == 1]
    return even, np.array([central(L) for L in even])


def _warm(L: EllipticL) -> None:
    central(L)
    hardy_z(1.0, L)


def family(
    E: Curve, lo: int, hi: int, tol: float = 1e-6, workers: int | None = None
) -> list[Twist]:
    even, value = even_values(E, lo, hi)
    vanish = [L for L, c in zip(even, value, strict=True) if abs(c) < tol]
    if workers == 1 or len(vanish) < 2:
        return [rank_two(L) for L in vanish]
    with ProcessPoolExecutor(workers, initializer=_warm, initargs=(vanish[-1],)) as pool:
        return list(pool.map(rank_two, vanish, chunksize=4))


def matrix_size(fam: list[Twist]) -> float:
    return float(np.mean([t.log_A for t in fam]))


def first_zeros(fam: list[Twist]) -> np.ndarray:
    return np.array([t.x[0] for t in fam])


def compare(
    x: np.ndarray, n: int, edges: np.ndarray, grid: range = range(2, 16)
) -> dict[str, tuple[float, float, int]]:
    out = {}
    for name, make in KERNELS.items():
        k: Kernel = make(n)
        best = fit_n(x, grid, make)[0]
        out[name] = (gap_chi2(x, k, edges), gap_ks(x, k), int(best))
    return out


def ks2(a: np.ndarray, b: np.ndarray) -> float:
    a, b = np.sort(a), np.sort(b)
    grid = np.concatenate([a, b])
    fa = np.searchsorted(a, grid, side="right") / len(a)
    fb = np.searchsorted(b, grid, side="right") / len(b)
    return float(np.max(np.abs(fa - fb)))


def excision_scan(
    x: np.ndarray, n: int, fractions: np.ndarray, size: int = 200_000, seed: int = 0
) -> tuple[np.ndarray, np.ndarray]:
    theta = jacobi_angles(n - 1, 1.5, size, seed)
    v = char_poly_at_one(theta)
    first = n * theta[:, 0] / math.pi
    cuts = np.quantile(v, fractions)
    return cuts, np.array([ks2(x, first[v >= c]) for c in cuts])
