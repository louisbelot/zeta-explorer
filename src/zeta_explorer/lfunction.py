from __future__ import annotations

import cmath
import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.characters import Character, factor
from zeta_explorer.zeros import illinois_many
from zeta_explorer.zeta import _B_OVER_FACT, _M_MAX, _log_sin, _out, loggamma, zeta

_CHUNK = 1 << 20


def dirichlet_l(s: ArrayLike, chi: Character, tol: float = 1e-16) -> complex | np.ndarray:
    z = np.array(s, dtype=complex, ndmin=1).ravel()
    if chi.is_principal:
        out = zeta(z, tol) * np.prod([1 - p ** (-z) for p, _ in factor(chi.q)], axis=0)
        return _out(np.reshape(out, np.shape(s)), s)
    prim = chi.primitive()
    out = _primitive_l(z, prim, tol)
    for p, _ in factor(chi.q):
        out *= 1 - prim(p) * p ** (-z)
    return _out(out.reshape(np.shape(s)), s)


def _primitive_l(s: np.ndarray, chi: Character, tol: float) -> np.ndarray:
    out = np.zeros_like(s)
    left = s.real < 0
    k = chi.parity
    trivial = left & (s.imag == 0) & ((s.real + k) % 2 == 0)
    right = ~left
    out[right] = _l_right(s[right], chi, tol)
    reflect = left & ~trivial
    if np.any(reflect):
        w = s[reflect]
        log_factor = (
            w * math.log(2 * math.pi) - math.log(math.pi) + (0.5 - w) * math.log(chi.q)
            + _log_sin(math.pi * (w + k) / 2) + loggamma(1 - w)
        )
        out[reflect] = chi.root_number() * np.exp(log_factor) * _l_right(1 - w, chi.conj(), tol)
    out[(s.real == 0) & (s.imag == 0) & (k == 0)] = 0
    if chi.is_real:
        out = np.where(s.imag == 0, out.real, out)
    return out


def _l_right(s: np.ndarray, chi: Character, tol: float) -> np.ndarray:
    q = chi.q
    a = np.flatnonzero(chi.num >= 0)
    c = chi.values[a]
    N = 10 + (np.abs(s.imag) / math.pi).astype(np.int64)

    out = _head(s, chi, N)
    w = N[:, None] + a[None, :] / q
    log_w = np.log(w)
    lead = (1 - s)[:, None] * log_w
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio = np.where(lead == 0, 1, np.expm1(lead) / lead)
    w_s = np.exp(-s[:, None] * log_w)
    tail = -ratio * log_w + w_s / 2

    poch = s[:, None].copy()
    power = w_s / w
    for j in range(1, _M_MAX + 1):
        term = _B_OVER_FACT[j] * poch * power
        tail += term
        if np.all(np.abs(term) < tol * 1e-2):
            break
        poch = poch * (s[:, None] + 2 * j - 1) * (s[:, None] + 2 * j)
        power = power / w**2
    return out + np.exp(-s * math.log(q)) * (tail @ c)


def _head(s: np.ndarray, chi: Character, N: np.ndarray) -> np.ndarray:
    out = np.zeros(len(s), dtype=complex)
    if len(s) == 0:
        return out
    n = np.arange(1, int(N.max()) * chi.q + 1)
    keep = chi.num[n % chi.q] >= 0
    n = n[keep]
    c = chi.values[n % chi.q]
    log_n = np.log(n)
    rows = max(1, _CHUNK // len(n))
    for i in range(0, len(s), rows):
        sl = slice(i, i + rows)
        terms = np.exp(-s[sl, None] * log_n) * c
        terms[n > (N[sl, None] * chi.q)] = 0
        out[sl] = terms.sum(axis=1)
    return out


def l_theta(t: ArrayLike, chi: Character) -> float | np.ndarray:
    prim = chi.primitive()
    ts = np.asarray(t, dtype=float)
    g = loggamma((0.5 + prim.parity + 1j * ts) / 2)
    eps = cmath.phase(prim.root_number())
    out = ts / 2 * math.log(prim.q / math.pi) + np.imag(g) - eps / 2
    return float(out) if np.ndim(out) == 0 else out


def l_z(t: ArrayLike, chi: Character) -> float | np.ndarray:
    prim = chi.primitive()
    ts = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.real(np.exp(1j * l_theta(ts, prim)) * dirichlet_l(0.5 + 1j * ts, prim))
    return float(out[0]) if np.ndim(t) == 0 else out


def mean_spacing(t: float, q: int) -> float:
    return 2 * math.pi / max(1.0, math.log(q * max(abs(t), 1.0) / (2 * math.pi)))


def l_zeros(
    chi: Character, t_lo: float, t_hi: float, per_gap: int = 8, rtol: float = 1e-13
) -> np.ndarray:
    prim = chi.primitive()
    step = mean_spacing(max(abs(t_lo), abs(t_hi)), prim.q) / per_gap
    ts = np.linspace(t_lo, t_hi, max(2, math.ceil((t_hi - t_lo) / step)) + 1)
    z = l_z(ts, prim)
    i = np.flatnonzero(z[:-1] * z[1:] < 0)
    exact = ts[z == 0]
    roots = illinois_many(lambda x: l_z(x, prim), ts[i], ts[i + 1], z[i], z[i + 1], rtol,
                          atol=1e-15)
    return np.sort(np.concatenate([roots, exact]))


def count_zeros_box(chi: Character, T: float, sigma: tuple[float, float] = (-0.5, 2.0)) -> int:
    prim = chi.primitive()
    if prim.is_principal:
        raise ValueError("use turing.verify_zeros for zeta")
    lo, hi = sigma
    corners = [lo - 1j * T, hi - 1j * T, hi + 1j * T, lo + 1j * T, lo - 1j * T]
    winding = 0.0
    for a, b in zip(corners[:-1], corners[1:], strict=True):
        winding += _arg_change(prim, a, b)
    n = winding / (2 * math.pi)
    if abs(n - round(n)) > 0.05:
        raise RuntimeError(f"winding number {n:.3f} is not close to an integer")
    trivial = len(range(-prim.parity, math.floor(lo), -2))
    return round(n) - trivial


def _arg_change(chi: Character, a: complex, b: complex, max_step: float = 0.5) -> float:
    u = np.linspace(0.0, 1.0, max(9, math.ceil(abs(b - a) / 0.05)))
    vals = dirichlet_l(a + (b - a) * u, chi)
    for _ in range(30):
        d = np.angle(vals[1:] / vals[:-1])
        bad = np.flatnonzero(np.abs(d) > max_step)
        if len(bad) == 0:
            return float(d.sum())
        mid = (u[bad] + u[bad + 1]) / 2
        u = np.insert(u, bad + 1, mid)
        vals = np.insert(vals, bad + 1, dirichlet_l(a + (b - a) * mid, chi))
    raise RuntimeError("argument did not settle, L may vanish on the contour")


@dataclass(frozen=True)
class LVerification:
    T: float
    on_line: int
    in_box: int
    zeros: np.ndarray

    @property
    def ok(self) -> bool:
        return self.on_line == self.in_box


def verify_zeros(chi: Character, T: float, max_per_gap: int = 256) -> LVerification:
    in_box = count_zeros_box(chi, T)
    per_gap = 8
    zeros = l_zeros(chi, -T, T, per_gap)
    while len(zeros) < in_box and per_gap < max_per_gap:
        per_gap *= 2
        zeros = l_zeros(chi, -T, T, per_gap)
    return LVerification(T, len(zeros), in_box, zeros)
