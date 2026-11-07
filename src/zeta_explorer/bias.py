from __future__ import annotations

import cmath
import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import Character, characters, factor
from zeta_explorer.explicit import _CHUNK, EULER_GAMMA, ei, mobius, primes
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

    def pdf(self, u: ArrayLike) -> np.ndarray:
        h, t, v = self._grid()
        u = np.asarray(u, dtype=float)
        w = np.full(len(t), h)
        w[0] = h / 2
        c = np.cos(np.multiply.outer(self.bias - u.ravel(), t))
        return (c @ (w * v) / math.pi).reshape(u.shape)

    def sample(self, n: int, rng: np.random.Generator) -> np.ndarray:
        out = self.bias + math.sqrt(self.tail / 2) * rng.standard_normal(n)
        amp = self.amp
        rows = max(1, _CHUNK // max(len(amp), 1))
        for i in range(0, n, rows):
            s = slice(i, min(n, i + rows))
            theta = rng.uniform(0, 2 * np.pi, (s.stop - s.start, len(amp)))
            out[s] += np.cos(theta) @ amp
        return out

    def explicit(self, x: ArrayLike) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        L = np.log(x).ravel()
        rho = 0.5 + 1j * self.gamma
        shift = 1j * np.pi * np.sign(self.gamma)
        out = np.zeros(len(L))
        rows = max(1, _CHUNK // max(len(rho), 1))
        for i in range(0, len(L), rows):
            s = slice(i, i + rows)
            out[s] = ((ei(L[s, None] * rho) - shift) @ self.coef).real
        for w, chi in self.chars:
            out += (w * _pi_main(L, chi)).real
        return (out * L / np.exp(L / 2)).reshape(x.shape)


def _log_at_zero(chi: Character) -> complex:
    k = chi.parity
    a = np.arange(1, chi.q)
    v = chi.values[1:]
    if k:
        z = complex(-np.sum(a * v) / chi.q)
    else:
        z = complex(np.sum(v * np.array([math.lgamma(t) for t in a / chi.q])))
    sigma = np.linspace(4.0, 0.01, 800)
    f = np.asarray(dirichlet_l(sigma + 0j, chi), dtype=complex) / sigma ** (1 - k)
    turns = round((np.unwrap(np.angle(f))[-1] - cmath.phase(z)) / (2 * math.pi))
    return complex(math.log(abs(z)), cmath.phase(z) + 2 * math.pi * turns)


def _prime_powers(chi: Character, k: int, y_max: float) -> tuple[np.ndarray, np.ndarray]:
    n, w = [], []
    for p in primes(int(y_max)).tolist():
        v = chi(p) ** k
        m, pm = 1, p
        while pm <= y_max:
            n.append(pm)
            w.append(v**m / m)
            m, pm = m + 1, pm * p
    order = np.argsort(n)
    return np.asarray(n, dtype=float)[order], np.cumsum(np.asarray(w, dtype=complex)[order])


def _pi_main(L: np.ndarray, chi: Character) -> np.ndarray:
    prim = chi.primitive()
    out = np.full(len(L), _log_at_zero(prim))
    if prim.parity == 0:
        out -= EULER_GAMMA + np.log(L)
    for p, _ in factor(chi.q):
        m = np.floor(L / math.log(p) + 1e-12).astype(np.int64)
        powers = np.cumsum([prim(p) ** j / j for j in range(1, m.max() + 1)])
        out -= np.where(m > 0, np.concatenate([[0], powers])[m], 0)
    k = 2
    while (live := L / k > math.log(2)).any():
        if mu := mobius(k):
            y = np.exp(L[live] / k)
            n, Pi = _prime_powers(chi, k, y.max())
            i = np.searchsorted(n, y * (1 + 1e-12), side="right") - 1
            out[live] += mu / k * np.where(i >= 0, Pi[np.maximum(i, 0)], 0)
        k += 1
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


def prime_race(q: int, a: int, b: int, n: int) -> tuple[np.ndarray, np.ndarray]:
    p = primes(n)
    d = np.cumsum((p % q == a % q).astype(np.int64) - (p % q == b % q))
    return p, d


def observed_density(p: np.ndarray, d: np.ndarray) -> float:
    span = np.diff(np.log(p))
    return float(np.sum(span[d[:-1] > 0]) / np.sum(span))
