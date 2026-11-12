from __future__ import annotations

import functools
import math
from collections.abc import Callable

import mpmath
import numpy as np
from numpy.typing import ArrayLike

Kernel = Callable[[np.ndarray, np.ndarray], np.ndarray]


def sine(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.sinc(x - y)


def symplectic(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.sinc(x - y) - np.sinc(x + y)


def even_orthogonal(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    return np.sinc(x - y) + np.sinc(x + y)


def _dirichlet(z: np.ndarray, m: float, n: float) -> np.ndarray:
    a = math.pi * z / (2 * n)
    den = np.sin(a)
    small = np.abs(den) < 1e-14
    out = np.sin(m * a) / (2 * n * np.where(small, 1.0, den))
    return np.where(small, m / (2 * n) * np.cos(m * a) / np.cos(a), out)


def symplectic_n(n: float) -> Kernel:
    if n < 1:
        raise ValueError("n must be >= 1")
    m = 2 * n + 1
    return lambda x, y: _dirichlet(x - y, m, n) - _dirichlet(x + y, m, n)


def even_orthogonal_n(n: float) -> Kernel:
    if n < 1:
        raise ValueError("n must be >= 1")
    m = 2 * n - 1
    return lambda x, y: _dirichlet(x - y, m, n) + _dirichlet(x + y, m, n)


def odd_orthogonal_n(n: float) -> Kernel:
    if n < 1:
        raise ValueError("n must be >= 1")
    m = 2 * n
    return lambda x, y: _dirichlet(x - y, m, n) - _dirichlet(x + y, m, n)


def odd_orthogonal_half(n: float) -> Kernel:
    k = odd_orthogonal_n(n)
    c = n / (n + 0.5)
    return lambda x, y: c * k(c * x, c * y)


def _jacobi_recurrence(n: int, alpha: float, beta: float) -> tuple[np.ndarray, np.ndarray, float]:
    ab = alpha + beta
    k = np.arange(n, dtype=float)
    b = np.empty(n)
    b[0] = (beta - alpha) / (ab + 2)
    b[1:] = (beta**2 - alpha**2) / ((2 * k[1:] + ab) * (2 * k[1:] + ab + 2))
    a = np.zeros(n)
    if n > 1:
        a[1] = 2 / (2 + ab) * math.sqrt((1 + alpha) * (1 + beta) / (3 + ab))
    j = k[2:]
    a[2:] = np.sqrt(4 * j * (j + alpha) * (j + beta) * (j + ab)
                    / ((2 * j + ab) ** 2 * (2 * j + ab + 1) * (2 * j + ab - 1)))
    log_mu = (ab + 1) * math.log(2) + math.lgamma(alpha + 1) + math.lgamma(beta + 1)
    log_mu -= math.lgamma(ab + 2)
    return b, a, math.exp(-log_mu / 2)


def jacobi_n(n: int, alpha: float, beta: float = -0.5, scale: float | None = None) -> Kernel:
    if n < 1 or n != int(n):
        raise ValueError("n must be a positive integer")
    if alpha <= -1 or beta <= -1:
        raise ValueError("the Jacobi weight needs alpha, beta > -1")
    b, a, p0 = _jacobi_recurrence(int(n), alpha, beta)
    c = math.pi / (n if scale is None else scale)

    def angle(x):
        h = c * np.asarray(x, dtype=float) / 2
        w = 2 ** (alpha + beta + 1) * np.sin(h) ** (2 * alpha + 1) * np.cos(h) ** (2 * beta + 1)
        return np.cos(2 * h), np.sqrt(c * w)

    def kernel(x, y):
        cx, wx = angle(x)
        cy, wy = angle(y)
        px0, px = 0.0, np.full_like(cx, p0)
        py0, py = 0.0, np.full_like(cy, p0)
        total = px * py
        for k in range(len(b) - 1):
            px0, px = px, ((cx - b[k]) * px - a[k] * px0) / a[k + 1]
            py0, py = py, ((cy - b[k]) * py - a[k] * py0) / a[k + 1]
            total = total + px * py
        return wx * wy * total

    return kernel


def _besselj(nu: float, z: np.ndarray) -> np.ndarray:
    f = np.frompyfunc(lambda v: float(mpmath.besselj(nu, v)), 1, 1)
    return np.asarray(f(z), dtype=float)


def bessel(alpha: float) -> Kernel:
    if alpha <= -0.5:
        raise ValueError("bessel needs alpha > -1/2, even_orthogonal is alpha = -1/2")

    def kernel(x, y):
        z = math.pi * np.asarray(x, dtype=float)
        w = math.pi * np.asarray(y, dtype=float)
        jz, jz1 = _besselj(alpha, z), _besselj(alpha + 1, z)
        jw, jw1 = _besselj(alpha, w), _besselj(alpha + 1, w)
        den = z * z - w * w
        same = np.abs(den) < 1e-12 * (1 + z * z)
        off = math.pi * np.sqrt(z * w) * (z * jz1 * jw - w * jz * jw1) / np.where(same, 1, den)
        diag = math.pi * z / 2 * (jz * jz - jz1 * _besselj(alpha - 1, z))
        diag = np.where(z == 0, 0.0, diag)
        return np.where(same, np.broadcast_to(diag, same.shape), off)

    return kernel


def independent_n(n: int) -> Kernel:
    return jacobi_n(n - 1, -0.5, scale=n)


def interaction_n(n: int) -> Kernel:
    return jacobi_n(n - 1, 1.5, scale=n)


def conditioned_n(n: int) -> Kernel:
    return jacobi_n(n - 1, 1.0, scale=n)


@functools.cache
def _legendre(m: int) -> tuple[np.ndarray, np.ndarray]:
    return np.polynomial.legendre.leggauss(m)


def _nystrom(
    kernel: Kernel, a: float, b: float, m: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u, w = _legendre(m)
    x = (b - a) / 2 * u + (a + b) / 2
    sw = np.sqrt(w * (b - a) / 2)
    return x, sw, np.eye(m) - sw[:, None] * kernel(x[:, None], x[None, :]) * sw[None, :]


def fredholm_det(kernel: Kernel, a: float, b: float, m: int = 40) -> float:
    return float(np.linalg.det(_nystrom(kernel, a, b, m)[2]))


def _gap_and_density(kernel: Kernel, s: float, m: int) -> tuple[float, float]:
    if s == 0:
        return 1.0, float(kernel(np.zeros(1), np.zeros(1))[0])
    x, sw, I_A = _nystrom(kernel, 0.0, s, m)
    edge = kernel(x, np.full(m, s))
    rho = np.linalg.solve(I_A, sw * edge)
    resolvent = float(kernel(np.array([s]), np.array([s]))[0]) + sw @ (edge * rho)
    E = float(np.linalg.det(I_A))
    return E, E * resolvent


def gap(kernel: Kernel, s: ArrayLike, m: int = 40) -> float | np.ndarray:
    ss = np.atleast_1d(np.asarray(s, dtype=float))
    out = np.array([_gap_and_density(kernel, si, m)[0] for si in ss])
    return float(out[0]) if np.ndim(s) == 0 else out


def first_density(kernel: Kernel, s: ArrayLike, m: int = 40) -> float | np.ndarray:
    ss = np.atleast_1d(np.asarray(s, dtype=float))
    out = np.array([_gap_and_density(kernel, si, m)[1] for si in ss])
    return float(out[0]) if np.ndim(s) == 0 else out


def first_density_mix(
    ns: ArrayLike, s: ArrayLike, family: Callable[[float], Kernel], groups: int = 12,
    m: int = 40,
) -> np.ndarray:
    ns = np.asarray(ns, dtype=float)
    ss = np.asarray(s, dtype=float)
    edges = np.linspace(ns.min(), ns.max(), groups + 1)
    idx = np.minimum(np.searchsorted(edges, ns, side="right") - 1, groups - 1)
    out = np.zeros_like(ss)
    for k in np.unique(idx):
        sel = idx == k
        n = float(ns[sel].mean())
        inside = ss < 0.98 * n
        out[inside] += sel.sum() * first_density(family(n), ss[inside], m)
    return out / len(ns)


def fit_n(
    first: ArrayLike, grid: ArrayLike, family: Callable[[float], Kernel] = symplectic_n,
    m: int = 40,
) -> tuple[float, float, float, np.ndarray]:
    x = np.asarray(first, dtype=float)
    ns = np.asarray(grid, dtype=float)
    s = np.linspace(0, x.max() * 1.001, 200)
    loglik = np.full(len(ns), -np.inf)
    for i, n in enumerate(ns):
        if n <= s[-1]:
            continue
        p = np.interp(x, s, first_density(family(n), s, m))
        if np.all(p > 0):
            loglik[i] = np.sum(np.log(p))
    best = int(np.argmax(loglik))
    inside = ns[loglik >= loglik[best] - 1.92]
    return float(ns[best]), float(inside.min()), float(inside.max()), loglik


def so_angles(
    n: int, size: int, seed: int = 0, batch: int = 10_000, odd: bool = False
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    dim = 2 * n + odd
    out = []
    for k in range(0, size, batch):
        Z = rng.standard_normal((min(batch, size - k), dim, dim))
        Q, R = np.linalg.qr(Z)
        Q *= np.sign(np.diagonal(R, axis1=1, axis2=2))[:, None, :]
        Q[np.linalg.det(Q) < 0, :, 0] *= -1
        theta = np.sort(np.abs(np.angle(np.linalg.eigvals(Q))), axis=1)
        out.append(theta[:, int(odd)::2])
    return np.concatenate(out)


def o_minus_angles(n: int, size: int, seed: int = 0, batch: int = 10_000) -> np.ndarray:
    rng = np.random.default_rng(seed)
    out = []
    for k in range(0, size, batch):
        Z = rng.standard_normal((min(batch, size - k), 2 * n, 2 * n))
        Q, R = np.linalg.qr(Z)
        Q *= np.sign(np.diagonal(R, axis1=1, axis2=2))[:, None, :]
        Q[np.linalg.det(Q) > 0, :, 0] *= -1
        theta = np.sort(np.abs(np.angle(np.linalg.eigvals(Q))), axis=1)
        out.append(theta[:, 1:-1:2])
    return np.concatenate(out)


def char_poly_at_one(theta: np.ndarray) -> np.ndarray:
    return np.prod(2 - 2 * np.cos(theta), axis=-1)


def excised_first(
    n: int, cutoff: float, size: int, seed: int = 0, odd: bool = False
) -> np.ndarray:
    theta = so_angles(n, size, seed, odd=odd)
    keep = char_poly_at_one(theta) >= cutoff
    return n * theta[keep, 0] / math.pi


def _jacobi_alpha(odd: bool) -> float:
    return 1.5 if odd else 0.5


def jacobi_moment(n: int, s: ArrayLike, alpha: float = -0.5) -> float | np.ndarray:
    a = alpha + 1
    ss = np.asarray(s, dtype=float)
    if np.any(ss <= -a):
        raise ValueError(f"the moment diverges for s <= {-a}")
    j = np.arange(n)[:, None]
    lg = np.vectorize(math.lgamma)
    log_m = n * ss.ravel() * math.log(4) + np.sum(
        lg(a + ss.ravel() + j) - lg(a + j) + lg(a + n + j - 0.5) - lg(a + n + j - 0.5 + ss.ravel()),
        axis=0,
    )
    out = np.exp(log_m).reshape(ss.shape)
    return float(out) if out.ndim == 0 else out


def so_moment(n: int, s: ArrayLike, odd: bool = False) -> float | np.ndarray:
    return jacobi_moment(n, s, 0.5 if odd else -0.5)


def small_value_law(n: int, odd: bool = False) -> tuple[float, float]:
    a = _jacobi_alpha(odd)
    log_r = -n * a * math.log(4) - math.lgamma(a)
    for j in range(n):
        log_r += math.lgamma(a + n + j - 0.5) - math.lgamma(n + j - 0.5)
        if j:
            log_r += math.lgamma(j) - math.lgamma(a + j)
    return math.exp(log_r), a


def binned_chi2(x: ArrayLike, model: ArrayLike, edges: ArrayLike) -> float:
    obs = np.histogram(x, edges)[0]
    expect = np.histogram(model, edges)[0] / len(model) * len(x)
    seen = expect > 0
    if np.any(obs[~seen] > 0):
        return math.inf
    return float(np.sum((obs[seen] - expect[seen]) ** 2 / expect[seen]))


def jacobi_angles(n: int, alpha: float, size: int, seed: int = 0, pool: int = 4) -> np.ndarray:
    if alpha == -0.5:
        return so_angles(n, size, seed)
    if alpha < 0.5:
        raise ValueError("reweighting SO(2n + 1) reaches alpha >= 1/2 only")
    theta = so_angles(n, pool * size, seed, odd=True)
    if alpha == 0.5:
        return theta[:size]
    w = char_poly_at_one(theta) ** (alpha - 0.5)
    pick = np.random.default_rng(seed + 1).choice(len(theta), size, p=w / w.sum())
    return theta[pick]


MODELS = {"independent": -0.5, "interaction": 1.5, "conditioned": 1.0}


def rank_two_first(
    n: int, model: str, size: int, seed: int = 0, cutoff: float = 0.0
) -> np.ndarray:
    theta = jacobi_angles(n - 1, MODELS[model], size, seed)
    keep = char_poly_at_one(theta) >= cutoff
    return n * theta[keep, 0] / math.pi


def gap_chi2(x: ArrayLike, kernel: Kernel, edges: ArrayLike, m: int = 40) -> float:
    obs = np.histogram(x, edges)[0]
    expect = -np.diff(gap(kernel, edges, m)) * len(np.asarray(x))
    return float(np.sum((obs - expect) ** 2 / expect))


def gap_ks(x: ArrayLike, kernel: Kernel, m: int = 40, grid: int = 300) -> float:
    xs = np.sort(np.asarray(x, dtype=float))
    s = np.linspace(0, xs[-1], grid)
    cdf = np.interp(xs, s, 1 - gap(kernel, s, m))
    k = np.arange(1, len(xs) + 1) / len(xs)
    return float(max(np.max(k - cdf), np.max(cdf - k + 1 / len(xs))))
