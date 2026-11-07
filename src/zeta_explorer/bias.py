from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import Character, characters
from zeta_explorer.explicit import _CHUNK, EULER_GAMMA
from zeta_explorer.lfunction import dirichlet_l, l_zeros


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


def _classes(q: int) -> list[Character]:
    out, seen = [], set()
    for chi in characters(q):
        if chi.is_principal:
            continue
        key = tuple(chi.values.round(12))
        if key in seen:
            continue
        seen.add(key)
        seen.add(tuple(chi.conj().values.round(12)))
        out.append(chi)
    return out


@dataclass(frozen=True)
class Race:
    q: int
    a: int
    b: int
    T: float
    bias: float
    gamma: np.ndarray
    coef: np.ndarray
    tail: float
    chars: tuple[tuple[complex, Character], ...]

    @property
    def amp(self) -> np.ndarray:
        return np.abs(self.coef) / np.abs(0.5 + 1j * self.gamma)

    @property
    def variance(self) -> float:
        return float(np.sum(self.amp**2) + self.tail) / 2

    def cf(self, t: ArrayLike) -> np.ndarray:
        t = np.asarray(t, dtype=float)
        return np.exp(-self.tail * t**2 / 4) * np.prod(bessel_j0(np.outer(self.amp, t)), axis=0)

    def _grid(self) -> tuple[float, np.ndarray, np.ndarray]:
        h = 2 * math.pi / (abs(self.bias) + 40 * math.sqrt(self.variance))
        t, vals = [np.zeros(1)], [np.ones(1)]
        k = 1
        while k < 200_000:
            block = h * np.arange(k, k + 512)
            v = self.cf(block)
            t.append(block)
            vals.append(v)
            k += 512
            if np.max(np.abs(v)) < 1e-17:
                break
        return h, np.concatenate(t), np.concatenate(vals)

    def density(self) -> float:
        h, t, v = self._grid()
        g = np.empty_like(t)
        g[0] = self.bias
        g[1:] = np.sin(self.bias * t[1:]) / t[1:] * v[1:]
        return 0.5 + h * (np.sum(g) - g[0] / 2) / math.pi

    def sample(self, n: int, rng: np.random.Generator) -> np.ndarray:
        out = self.bias + math.sqrt(self.tail / 2) * rng.standard_normal(n)
        amp = self.amp
        rows = max(1, _CHUNK // max(len(amp), 1))
        for i in range(0, n, rows):
            s = slice(i, min(n, i + rows))
            theta = rng.uniform(0, 2 * np.pi, (s.stop - s.start, len(amp)))
            out[s] += np.cos(theta) @ amp
        return out




def race(q: int, a: int, b: int, T: float = 1000.0) -> Race:
    if math.gcd(a, q) != 1 or math.gcd(b, q) != 1 or (a - b) % q == 0:
        raise ValueError("a and b must be distinct residues coprime to q")
    gamma, coef, tail, chars = [], [], 0.0, []
    for chi in _classes(q):
        delta = chi(a).conjugate() - chi(b).conjugate()
        weight = 1 if chi.is_real else 2
        if abs(delta) < 1e-12:
            continue
        chars.append((weight * delta, chi))
        lo = 0.0 if chi.is_real else -T
        g = l_zeros(chi, lo, T)
        gamma.append(g)
        coef.append(np.full(len(g), -2 * delta))
        full = zero_sum(chi) * (0.5 if chi.is_real else 1.0)
        tail += 4 * abs(delta) ** 2 * (full - np.sum(1 / (0.25 + g**2)))
    bias = float(square_roots(q, b) - square_roots(q, a))
    return Race(q, a, b, T, bias, np.concatenate(gamma), np.concatenate(coef), tail, tuple(chars))