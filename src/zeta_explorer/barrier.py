from __future__ import annotations

import math

import mpmath
import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.effective import _log_mt_mp, alpha, n_terms
from zeta_explorer.multieval import log_phases

X0 = 60_000_083_952.0
PRIMES = (2, 3, 5)
_CHUNK = 1 << 22
_GL = [(0.5 * (1 + u), 0.5 * w) for u, w in zip(*np.polynomial.legendre.leggauss(8), strict=True)]


def _ranges(lo: int, hi: int):
    for a in range(lo, hi + 1, _CHUNK):
        yield np.arange(a, min(hi, a + _CHUNK - 1) + 1, dtype=np.int64)


def log_b(n: ArrayLike, t: float) -> np.ndarray:
    return t / 4 * np.log(n) ** 2


def bsum(M: int, sigma: float, t: float) -> float:
    return sum(float(np.exp(log_b(n, t) - sigma * np.log(n)).sum()) for n in _ranges(1, M))


def lemma_tail(M: int, N: float, sigma: float, t: float) -> float:
    if sigma <= t / 2 * math.log(N):
        raise ValueError("Lemma 8.2 needs sigma > (t/2) log N")
    ends = [(1 - sigma) * math.log(v) + t / 4 * math.log(v) ** 2 for v in (M, N)]
    return math.exp(max(ends)) * math.log(N / M)


def erfi_tail(M: int, N: float, sigma: float, t: float) -> float:
    with mpmath.workdps(30):
        r = mpmath.sqrt(t)

        def g(v):
            return mpmath.erfi(r / 2 * mpmath.log(v) - (sigma - 1) / r)

        return float(mpmath.sqrt(mpmath.pi / t) * mpmath.exp(-((sigma - 1) ** 2) / t)
                     * (g(N) - g(M)))


def x_of(N: float, t: float) -> float:
    return 4 * math.pi * (N * N - t / 16)


def sigma_floor(N: float, y: float, t: float) -> float:
    xN = x_of(N, t)
    return (1 + y) / 2 + t / 4 * math.log(xN / (4 * math.pi)) - t / (2 * xN * xN)


def gamma_constant(N: float, y: float, t: float, n_max: float) -> float:
    return (math.exp(0.02 * y) * (1 - t / (16 * N * N)) ** (-y / 2)
            * n_max ** (t / (2 * (x_of(N, t) - 6))))


def far_bound(N: float, t: float, y: float, M: int) -> tuple[float, float]:
    a = sigma_floor(N, y, t)
    c = gamma_constant(N, y, t, N) * N**-y
    A = bsum(M, a, t) + c * bsum(M, a - y, t)
    B = lemma_tail(M, N, a, t) + c * lemma_tail(M, N, a - y, t)
    return A, B


def far_sup(N_lo: float, N_hi: float, t: float, y: float, M: int, samples: int = 400) -> float:
    return max(sum(far_bound(N, t, y, M)) for N in np.geomspace(N_lo, N_hi, samples))


def top_bound(N0: int, N1: int, t: float, y: float = 1.0) -> float:
    a = sigma_floor(N0, y, t)
    c = gamma_constant(N0, y, t, N1)
    head = bsum(N1, a, t)
    return head - 1 + c * (bsum(N0, a - y, t) * N0**-y + head - bsum(N0, a, t))


def mollifier_terms(t: float) -> list[tuple[int, float]]:
    out = []
    for mask in range(1 << len(PRIMES)):
        ps = [p for i, p in enumerate(PRIMES) if mask >> i & 1]
        out.append((math.prod(ps), math.prod(-math.exp(log_b(p, t)) for p in ps)))
    return sorted(out)


def mollifier_range(sigma: float, t: float) -> tuple[float, float, float]:
    r = [math.exp(log_b(p, t)) / p**sigma for p in PRIMES]
    return (math.prod(1 - v for v in r), math.prod(1 + v for v in r),
            sum(math.asin(v) for v in r))


def coefficients(n: np.ndarray, cap: ArrayLike, t: float,
                 y: float) -> tuple[np.ndarray, np.ndarray]:
    beta = np.zeros(len(n))
    lam = np.zeros(len(n))
    cap = np.broadcast_to(cap, n.shape)
    for d, a in mollifier_terms(t):
        m = (n % d == 0) & (n <= d * cap)
        q = (n[m] // d).astype(float)
        bq = a * np.exp(log_b(q, t))
        beta[m] += bq
        lam[m] += bq * q**y
    return beta, lam


def lower_edge(N_lo: int, N_hi: int, t: float = 0.2, y: float = 0.2, improved: bool = True,
               every_n: bool = True, c: float | None = None) -> float:
    c = gamma_constant(N_lo, y, t, N_hi) if c is None else c
    r = c * N_lo**-y
    k = (1 - r) / (1 + r)
    sigma = sigma_floor(N_lo, y, t)
    D = math.prod(PRIMES)
    divisors = [d for d, _ in mollifier_terms(t)]
    total = 0.0
    for n in _ranges(2, D * N_hi):
        caps = [float(N_lo), *(n / d for d in divisors)] if every_n else [float(N_hi)]
        worst = np.zeros(len(n))
        for cap in caps:
            beta, lam = coefficients(n, cap, t, y)
            if improved:
                term = np.maximum(np.abs(beta - r * lam), k * np.abs(beta + r * lam))
            else:
                term = np.abs(beta) + r * np.abs(lam)
            if every_n:
                term = np.where((cap >= N_lo) & (cap <= N_hi), term, 0.0)
            worst = np.maximum(worst, term)
        total += float((worst * np.exp(-sigma * np.log(n))).sum())
    return 1 - r - total


def polar_distance(mag: float, arg: float) -> float:
    if arg >= math.pi:
        return 0.0
    return mag if arg <= math.pi / 2 else mag * math.sin(arg)


def _tail(z: float, k: int) -> float:
    term = z**k / math.factorial(k)
    total = 0.0
    while term > 1e-300:
        total += term
        k += 1
        term *= z / k
        if term < 1e-30 * total:
            break
    return total


class Barrier:
    def __init__(self, X: float = X0, y_lo: float = 0.2, y_hi: float = 1.0, width: float = 1.0,
                 t_hi: float = 0.2, block: float = 0.01, tol: float = 1e-13):
        Ns = {int(n_terms(X + e, t)) for e in (-width / 2, width / 2) for t in (0.0, t_hi)}
        if len(Ns) != 1:
            raise ValueError("N must be constant on the barrier")
        self.X, self.y_lo, self.y_hi, self.width, self.t_hi = X, y_lo, y_hi, width, t_hi
        self.N = Ns.pop()
        self.block, self.tol = block, tol
        self.yc = (y_lo + y_hi) / 2
        self.L = np.log(np.arange(1, self.N + 1, dtype=float))
        self.P = np.exp(1j * log_phases(X / 2, self.N))
        with mpmath.workdps(30):
            s = (1 + mpmath.mpf(self.yc) - 1j * mpmath.mpf(X)) / 2
            d = _log_mt_mp(1 - s, 0) - _log_mt_mp(s, 0)
            self.d0 = complex(float(d.real), float(d.imag % (2 * mpmath.pi)))
        self._blocks: dict[int, tuple] = {}

    def _s(self, u: np.ndarray, y: np.ndarray) -> np.ndarray:
        return (1 + y - 1j * (self.X + u)) / 2

    def log_gamma(self, u: ArrayLike, y: ArrayLike, t: float) -> np.ndarray:
        u, y = np.asarray(u, dtype=float), np.asarray(y, dtype=float)
        dz = u + 1j * (y - self.yc)
        acc = np.zeros(np.broadcast(u, y).shape, dtype=complex)
        for xi, w in _GL:
            s = self._s(xi * u, self.yc + xi * (y - self.yc))
            acc += w * (alpha(1 - s) + alpha(s))
        s = self._s(u, y)
        return self.d0 + 0.5j * dz * acc + t / 4 * (alpha(1 - s) ** 2 - alpha(s) ** 2)

    def omegas(self, u: ArrayLike, y: ArrayLike, t: float) -> tuple[np.ndarray, np.ndarray]:
        u, y = np.asarray(u, dtype=float), np.asarray(y, dtype=float)
        s = self._s(u, y)
        return (1 + y - 1j * u) / 2 + t / 2 * alpha(s), (1 - y + 1j * u) / 2 + t / 2 * alpha(1 - s)

    def _orders(self, lam: float, mu: float, weight: float) -> tuple[int, int]:
        J = 1
        while _tail(lam, J) * math.exp(mu) * weight > self.tol / 2:
            J += 1
        K = 1
        while math.exp(lam) * _tail(mu, K) * weight > self.tol / 2:
            K += 1
        return J, K

    def corners(self, t1: float, t2: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        g = np.array([(u, y, t) for u in (-self.width / 2, self.width / 2)
                      for y in (self.y_lo, self.y_hi) for t in (t1, t2)])
        return g[:, 0], g[:, 1], g[:, 2]

    def _span(self, t: float) -> tuple[int, float]:
        count = math.ceil(self.t_hi / self.block - 1e-9)
        i = min(math.floor(t / self.block + 1e-9), count - 1)
        return i, min((i + 1) * self.block, self.t_hi)

    def _block(self, i: int):
        if i in self._blocks:
            return self._blocks[i]
        t1, t2 = i * self.block, min((i + 1) * self.block, self.t_hi)
        tc = (t1 + t2) / 2
        w1c, w2c = (complex(w) for w in self.omegas(0.0, self.yc, tc))
        u, y, t = self.corners(t1, t2)
        s = self._s(u, y)
        g = float(np.exp(self.log_gamma(u, y, t).real.max()))
        amax = float(max(np.abs(alpha(s)).max(), np.abs(alpha(1 - s)).max()))
        Lmax = self.L[-1]
        hy = max(self.yc - self.y_lo, self.y_hi - self.yc)
        mu = (0.5 * math.hypot(hy, self.width / 2) + (t2 - t1) / 4 * amax + 1e-9) * Lmax
        lam = (t2 - t1) / 8 * Lmax * Lmax
        base = tc / 4 * self.L**2
        v1 = np.exp(base - w1c * self.L) * self.P
        v2 = np.exp(base - w2c * self.L) * np.conj(self.P)
        J, K = self._orders(lam, mu, float(np.abs(v1).sum() + g * np.abs(v2).sum()))
        S1 = np.empty(2 * J + K, dtype=complex)
        S2 = np.empty(2 * J + K, dtype=complex)
        for m in range(2 * J + K):
            S1[m], S2[m] = v1.sum(), v2.sum()
            v1 = v1 * self.L
            v2 = v2 * self.L
        self._blocks[i] = (tc, w1c, w2c, J, K, S1, S2)
        return self._blocks[i]

    @staticmethod
    def _series(S: np.ndarray, J: int, K: int, a: float, b: np.ndarray) -> np.ndarray:
        B = np.empty((K, len(b)), dtype=complex)
        B[0] = 1
        for k in range(1, K):
            B[k] = B[k - 1] * b / k
        out = np.zeros(len(b), dtype=complex)
        aj = 1.0
        for j in range(J):
            out += aj * (S[2 * j : 2 * j + K] @ B)
            aj *= a / (j + 1)
        return out

    def f(self, u: ArrayLike, y: ArrayLike, t: float) -> np.ndarray:
        u, y = np.broadcast_arrays(np.atleast_1d(np.asarray(u, dtype=float)),
                                   np.atleast_1d(np.asarray(y, dtype=float)))
        tc, w1c, w2c, J, K, S1, S2 = self._block(self._span(t)[0])
        w1, w2 = self.omegas(u, y, t)
        a = (t - tc) / 4
        g = np.exp(self.log_gamma(u, y, t))
        return self._series(S1, J, K, a, w1c - w1) + g * self._series(S2, J, K, a, w2c - w2)
