from __future__ import annotations

import math

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import Character
from zeta_explorer.explicit import _CHUNK, EULER_GAMMA
from zeta_explorer.lfunction import dirichlet_l


def bessel_j0(z: ArrayLike) -> np.ndarray:
    z = np.abs(np.asarray(z, dtype=float))
    top = float(z.max(initial=0.0))
    m = math.ceil((top + 10 * top ** (1 / 3) + 30) / 2)
    theta = np.pi * (np.arange(m) + 0.5) / m
    flat = z.ravel()
    out = np.empty_like(flat)
    rows = max(1, _CHUNK // m)
    for i in range(0, len(flat), rows):
        s = slice(i, i + rows)
        out[s] = np.cos(flat[s, None] * np.sin(theta)).mean(axis=1)
    return out.reshape(z.shape)


def square_roots(q: int, a: int) -> int:
    x = np.arange(q, dtype=np.int64)
    return int(np.count_nonzero(x * x % q == a % q))


def l_log_derivative_one(chi: Character, radius: float = 0.4, nodes: int = 48) -> complex:
    if chi.is_principal:
        raise ValueError("L(s, χ0) has a pole at s = 1")
    w = np.exp(2j * np.pi * np.arange(nodes) / nodes)
    vals = dirichlet_l(1 + radius * w, chi)
    L1 = complex(dirichlet_l(1.0 + 0j, chi))
    dL = complex(np.mean(vals / w)) / radius
    return dL / L1


def zero_sum(chi: Character) -> float:
    prim = chi.primitive()
    digamma = -EULER_GAMMA - (2 * math.log(2) if prim.parity == 0 else 0.0)
    lam = math.log(prim.q / math.pi) + digamma
    return lam + 2 * l_log_derivative_one(prim).real