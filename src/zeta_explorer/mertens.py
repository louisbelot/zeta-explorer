from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from flint import acb, acb_series, arb
from numpy.typing import ArrayLike

from zeta_explorer.certified import Certificate, certify, precision, z_sign
from zeta_explorer.riemann_siegel import z_rs
from zeta_explorer.zeros import illinois


def mobius_sieve(n: int) -> np.ndarray:
    mu = np.ones(n + 1, dtype=np.int8)
    mu[0] = 0
    is_prime = np.ones(n + 1, dtype=bool)
    is_prime[:2] = False
    for p in range(2, math.isqrt(n) + 1):
        if is_prime[p]:
            is_prime[p * p :: p] = False
    for p in np.flatnonzero(is_prime):
        mu[p::p] *= -1
    for p in np.flatnonzero(is_prime[: math.isqrt(n) + 1]):
        mu[p * p :: p * p] = 0
    return mu


def mertens_table(n: int) -> np.ndarray:
    return np.cumsum(mobius_sieve(n), dtype=np.int64)


def mertens(x: int, u: int | None = None) -> int:
    if x < 1:
        return 0
    if u is None:
        u = max(math.isqrt(x) + 1, int(x ** (2 / 3) / 4))
    small = mertens_table(max(u, math.isqrt(x) + 1))
    if x < len(small):
        return int(small[x])
    K = x // (len(small))
    L = np.zeros(K + 1, dtype=np.int64)
    for k in range(K, 0, -1):
        v = x // k
        r = math.isqrt(v)
        s = v // (r + 1)
        total = 1
        nb = min(r, K // k)
        if nb >= 2:
            total -= int(L[k * np.arange(2, nb + 1)].sum())
        d = np.arange(max(nb, 1) + 1, r + 1, dtype=np.int64)
        total -= int(small[v // d].sum())
        q = np.arange(1, s + 1, dtype=np.int64)
        counts = v // q - v // (q + 1)
        counts[-1] = v // s - r
        total -= int((counts * small[q]).sum())
        L[k] = total
    return int(L[1])


@dataclass(frozen=True)
class Zero:
    gamma: arb
    deriv: acb
    a: arb
    psi: arb


def zeta_deriv(s: acb) -> tuple[acb, acb]:
    f = acb_series([s, 1], prec=2).zeta()
    return f.coeffs()[0], f.coeffs()[1]


def refine(t: float, bits: int) -> arb:
    g = arb(t)
    b = 32
    while b < bits:
        b = min(2 * b, bits)
        with precision(b + 32):
            f, d = zeta_deriv(acb(0.5, g))
            g = (g + (-f / d).imag).mid()
    with precision(bits + 32):
        f, d = zeta_deriv(acb(0.5, g))
        return (g + (-f / d).imag).mid()


def isolate(t: float, bits: int) -> Zero:
    with precision(bits + 32):
        g = refine(t, bits)
        d = arb(2) ** (-bits) * (1 + abs(g))
        top = 4 * (bits + 32)
        lo, hi = z_sign(g - d, bits + 32, top), z_sign(g + d, bits + 32, top)
        if lo * hi != -1:
            raise RuntimeError(f"no certified sign change around t = {t}")
        return at(arb(g, d), bits)


def at(gamma: arb, bits: int) -> Zero:
    with precision(bits + 32):
        rho = acb(0.5, gamma)
        deriv = zeta_deriv(rho)[1]
        w = rho * deriv
        return Zero(gamma, deriv, 1 / abs(w), w.arg())


def brackets(cert: Certificate) -> list[tuple[float, float]]:
    out = []
    for b in cert.blocks:
        for i in range(len(b.signs) - 1):
            if b.signs[i] * b.signs[i + 1] < 0:
                out.append((b.points[i], b.points[i + 1]))
    return out


def certified_zeros(T: float, bits: int = 256) -> list[Zero]:
    cert = certify(T)
    out = []
    for lo, hi in brackets(cert):
        if lo > T:
            break
        t = illinois(z_rs, lo, hi, float(z_rs(lo)), float(z_rs(hi)), tol=1e-11)
        z = isolate(t, bits)
        if not (lo < z.gamma < hi):
            raise RuntimeError(f"refined zero left its bracket ({lo}, {hi})")
        out.append(z)
    return [z for z in out if z.gamma < T]


def kernel(u: arb) -> arb:
    return (1 - u) * (arb.pi() * u).cos() + (arb.pi() * u).sin() / arb.pi()


def kernel_np(u: ArrayLike) -> np.ndarray:
    u = np.asarray(u, dtype=float)
    return (1 - u) * np.cos(np.pi * u) + np.sin(np.pi * u) / np.pi


def work_bits(y: arb) -> int:
    return int(float((abs(y) + 1).log()) / math.log(2)) + 128


def ingham(y: arb, zeros: list[Zero]) -> arb:
    with precision(work_bits(y)):
        top = zeros[-1].gamma
        total = arb(0)
        for z in zeros:
            total += z.a * kernel(z.gamma / top) * (z.gamma * y - z.psi).cos()
        return 2 * total


def arrays(zeros: list[Zero]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    g = np.array([float(z.gamma.mid()) for z in zeros])
    a = np.array([float(z.a.mid()) for z in zeros])
    psi = np.array([float(z.psi.mid()) for z in zeros])
    return g, a, psi


def trivial(x: ArrayLike, terms: int = 30) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.full_like(x, -2.0)
    for n in range(1, terms + 1):
        z = float(acb(2 * n + 1).zeta().real)
        out += (-1) ** (n - 1) * np.exp(2 * n * np.log(2 * np.pi / x) - math.lgamma(2 * n + 1)) / (
            n * z
        )
    return out


def q_from_zeros(x: ArrayLike, zeros: list[Zero], smooth: bool = False) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    g, a, psi = arrays(zeros)
    w = kernel_np(g / g[-1]) if smooth else np.ones_like(g)
    y = np.log(x)
    q = 2 * (np.cos(np.outer(y, g) - psi) @ (a * w))
    return q + trivial(x) / np.sqrt(x)
