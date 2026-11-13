from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import ArrayLike

from zeta_explorer.afe import afe_length, l_afe
from zeta_explorer.bias import Race, prime_race, race
from zeta_explorer.characters import conrey, primitive_characters, quadratic, totient
from zeta_explorer.continuation import dirichlet_sums, em_sums, eta_sums
from zeta_explorer.delta import SatoTate
from zeta_explorer.elliptic import Curve, EllipticL, hardy_z, zeros_below
from zeta_explorer.explicit import chebyshev_psi, pi_from_zeros, prime_pi, psi_from_zeros, riemann_r
from zeta_explorer.gram import good_gram, z_values
from zeta_explorer.heat import ht
from zeta_explorer.lfunction import dirichlet_l, l_z, l_zeros
from zeta_explorer.li import li_trend
from zeta_explorer.lowlying import OneLevel, w_symplectic, w_unitary
from zeta_explorer.moments import keating_snaith, polynomial_integral
from zeta_explorer.multieval import ZWindow
from zeta_explorer.ratios import quadratic_density
from zeta_explorer.riemann_siegel import theta, z_exact, z_rs
from zeta_explorer.rmt import (
    bessel,
    binned_chi2,
    char_poly_at_one,
    conditioned_n,
    even_orthogonal,
    even_orthogonal_n,
    first_density,
    first_density_mix,
    fit_n,
    gap,
    independent_n,
    interaction_n,
    odd_orthogonal_half,
    odd_orthogonal_n,
    symplectic,
    symplectic_n,
)
from zeta_explorer.spacing import (
    first_zeros,
    gue_spacing,
    montgomery,
    normalized_spacings,
    pair_correlation,
    unfold,
    wigner_surmise,
)
from zeta_explorer.values import (
    arithmetic_cumulants,
    cue_cumulants,
    cue_density,
    hybrid_density,
    sample_cumulants,
)
from zeta_explorer.zeros import find_zeros, gram_points
from zeta_explorer.zeta import zeta

if TYPE_CHECKING:
    from matplotlib.figure import Figure

    from zeta_explorer.symfamily import Family as SymFamily

SERIES = "#2a78d6"
ACCENT = "#eb6834"
ORDINAL = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"


def plot_z(t_min: float = 0.0, t_max: float = 60.0, samples: int = 1500) -> Figure:
    import matplotlib.pyplot as plt

    ts = np.linspace(t_min, t_max, samples)
    zs = [z_exact(t) for t in ts]
    zeros = find_zeros(max(t_min, 10.0), t_max)

    fig, ax = plt.subplots(figsize=(10, 4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    ax.plot(ts, zs, color=SERIES, linewidth=2)
    ax.plot(zeros, [0] * len(zeros), "o", color=SERIES, markersize=8,
            markeredgecolor=SURFACE, markeredgewidth=2)

    _style(ax)
    ax.set_xlim(t_min, t_max)
    ax.set_xlabel("t", color=INK)
    ax.set_ylabel("Z(t)", color=INK)
    ax.set_title(
        f"Riemann-Siegel Z(t): every sign change is a zero of ζ(1/2 + it)  "
        f"({len(zeros)} zeros up to t = {t_max:g})",
        color=INK, fontsize=12, loc="left",
    )
    fig.tight_layout()
    return fig


def plot_gram(t_min: float = 262.0, t_max: float = 300.0, samples: int = 1500) -> Figure:
    import matplotlib.pyplot as plt

    ts = np.linspace(t_min, t_max, samples)
    n_lo = int(np.ceil(theta(t_min) / np.pi))
    n_hi = int(np.floor(theta(t_max) / np.pi))
    ns = np.arange(n_lo, n_hi + 1)
    g = gram_points(ns)
    z = z_values(g)
    good = good_gram(z, n0=n_lo)

    fig, ax = plt.subplots(figsize=(10, 4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    ax.plot(ts, z_values(ts), color=SERIES, linewidth=2)
    for mask, color, label in ((good, SERIES, "good"), (~good, ACCENT, "bad")):
        ax.plot(g[mask], z[mask], "o", color=color, markersize=8, markeredgecolor=SURFACE,
                markeredgewidth=2, label=f"{label} Gram point", linestyle="none")
    for n, t, zt in zip(ns[~good], g[~good], z[~good], strict=True):
        ax.annotate(f"g_{n}", (t, zt), xytext=(8, -14), textcoords="offset points",
                    color=INK_SECONDARY, fontsize=9)

    _style(ax)
    ax.set_xlim(t_min, t_max)
    ax.set_xlabel("t", color=INK)
    ax.set_ylabel("Z(t)", color=INK)
    ax.legend(frameon=False, fontsize=9, loc="lower left", labelcolor=INK_SECONDARY)
    ax.set_title("Gram points: g_n is good when (-1)^n Z(g_n) > 0", color=INK, fontsize=12,
                 loc="left")
    fig.tight_layout()
    return fig


def plot_rs_error(t_min: float = 200.0, t_max: float = 20_000.0, bins: int = 40) -> Figure:
    import matplotlib.pyplot as plt

    edges = np.geomspace(t_min, t_max, bins + 1)
    ts = np.concatenate(
        [np.linspace(a, b, 25, endpoint=False) for a, b in zip(edges[:-1], edges[1:], strict=True)]
    )
    exact = np.array([z_exact(t) for t in ts])
    centers = np.sqrt(edges[:-1] * edges[1:])

    fig, ax = plt.subplots(figsize=(10, 5), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    for k, color in enumerate(ORDINAL):
        err = np.abs(z_rs(ts, corrections=k) - exact).reshape(bins, 25).max(axis=1)
        ax.loglog(centers, err, color=color, linewidth=2)
        label = "C0" if k == 0 else f"C0 to C{k}"
        ax.annotate(f"{label}, $t^{{-{(2 * k + 3) / 4:g}}}$", (centers[0], err[0]),
                    xytext=(-8, 0), textcoords="offset points", ha="right", va="center",
                    color=INK_SECONDARY, fontsize=9)

    _style(ax)
    ax.set_xlim(t_min / 3.5, t_max)
    ax.set_xlabel("t", color=INK)
    ax.set_ylabel("max |error| in bin", color=INK)
    ax.set_title(
        "Riemann-Siegel error with more correction terms (floor: double-precision phases)",
        color=INK, fontsize=12, loc="left",
    )
    fig.tight_layout()
    return fig


def phase_colors(w: np.ndarray) -> np.ndarray:
    from matplotlib.colors import hsv_to_rgb

    hue = (np.angle(w) / (2 * np.pi)) % 1
    with np.errstate(divide="ignore", invalid="ignore"):
        bands = np.log2(np.abs(w)) % 1
    value = np.where(np.isfinite(bands), 0.62 + 0.38 * bands, 1.0)
    return hsv_to_rgb(np.stack([hue, np.full_like(hue, 0.78), value], axis=-1))


def plot_domain(
    t_min: float = 0.0, t_max: float = 60.0, sigma_min: float = -5.0, sigma_max: float = 6.0,
    width: int = 1800,
) -> Figure:
    import matplotlib.pyplot as plt

    height = int(width * (sigma_max - sigma_min) / (t_max - t_min))
    t = t_min + (np.arange(width) + 0.5) * (t_max - t_min) / width
    sigma = sigma_min + (np.arange(height) + 0.5) * (sigma_max - sigma_min) / height
    s = sigma[:, None] + 1j * t[None, :]
    w = np.full(s.shape, np.inf + 0j)
    finite = s != 1
    w[finite] = zeta(s[finite])

    fig, ax = plt.subplots(figsize=(12, 12 * height / width + 1.2), facecolor=SURFACE)
    ax.imshow(phase_colors(w), origin="lower", extent=(t_min, t_max, sigma_min, sigma_max),
              aspect="auto", interpolation="lanczos")
    for line in (0, 1):
        ax.axhline(line, color=SURFACE, linewidth=0.8, alpha=0.7)
    ax.axhline(0.5, color=SURFACE, linewidth=1, linestyle=(0, (4, 4)))
    for side in ax.spines.values():
        side.set_visible(False)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.set_xlabel("t = Im s", color=INK)
    ax.set_ylabel("σ = Re s", color=INK)
    ax.set_title(
        "Phase portrait of ζ(σ + it): hue is arg ζ, shading steps at every doubling of |ζ|",
        color=INK, fontsize=12, loc="left",
    )
    fig.tight_layout()
    return fig


def plot_spiral(t_max: float = 40.0, samples: int = 4000) -> Figure:
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    t = np.linspace(0, t_max, samples)
    w = zeta(0.5 + 1j * t)
    points = np.stack([w.real, w.imag], axis=1)
    segments = np.stack([points[:-1], points[1:]], axis=1)
    zeros = find_zeros(10, t_max)

    fig, ax = plt.subplots(figsize=(8.5, 6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    lines = LineCollection(segments, cmap="viridis", linewidth=1.6)
    lines.set_array((t[:-1] + t[1:]) / 2)
    ax.add_collection(lines)
    ax.plot([w[0].real], [w[0].imag], "o", color=INK, markersize=6)
    ax.annotate("t = 0, ζ(1/2) = −1.4604", (w[0].real, w[0].imag), xytext=(8, -14),
                textcoords="offset points", color=INK_SECONDARY, fontsize=9)
    ax.plot([0], [0], "o", color=ACCENT, markersize=7, markeredgecolor=SURFACE,
            markeredgewidth=1.5)
    ax.annotate(f"origin: the curve passes here {len(zeros)} times", (0, 0), xytext=(-10, 16),
                textcoords="offset points", ha="right", color=ACCENT, fontsize=9)
    _style(ax)
    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.set_xlabel("Re ζ(1/2 + it)", color=INK)
    ax.set_ylabel("Im ζ(1/2 + it)", color=INK)
    bar = fig.colorbar(lines, ax=ax, shrink=0.7, pad=0.02)
    bar.set_label("t", color=INK)
    bar.outline.set_visible(False)
    bar.ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.set_title(f"ζ(1/2 + it) for 0 ≤ t ≤ {t_max:g}", color=INK, fontsize=12, loc="left")
    fig.tight_layout()
    return fig


def plot_continuation(s: complex = 0.5 + 14.134725141734693j, N: int = 400) -> Figure:
    import matplotlib.pyplot as plt

    target = zeta(s)
    series = (
        ("Σ n^(−s), n ≤ N", dirichlet_sums(s, N), ORDINAL[0]),
        ("alternating η series / (1 − 2^(1−s))", eta_sums(s, N), ORDINAL[2]),
        ("Σ n^(−s) + Euler-Maclaurin tail", em_sums(s, N), ACCENT),
    )

    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 5), facecolor=SURFACE,
                                      gridspec_kw={"width_ratios": [1, 1.25]})
    n = np.arange(1, N + 1)
    for label, sums, color in series:
        left.plot(sums.real, sums.imag, color=color, linewidth=1.2, label=label)
        right.loglog(n, np.abs(sums - target), color=color, linewidth=1.6, label=label)
    left.plot([target.real], [target.imag], "o", color=INK, markersize=6)
    left.annotate("ζ(s)", (target.real, target.imag), xytext=(8, 6),
                  textcoords="offset points", color=INK, fontsize=10,
                  bbox={"facecolor": SURFACE, "edgecolor": "none", "pad": 1})
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    left.set_aspect("equal", adjustable="datalim")
    left.set_xlabel("real part", color=INK)
    left.set_ylabel("imaginary part", color=INK)
    left.set_title("Partial sums in the complex plane", color=INK, fontsize=11, loc="left")
    right.set_xlabel("N", color=INK)
    right.set_ylabel("|partial sum − ζ(s)|", color=INK)
    right.set_title("Distance to ζ(s)", color=INK, fontsize=11, loc="left")
    right.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="lower left")
    fig.suptitle(f"Three ways to sum ζ(s) at s = {s.real:g} + {s.imag:.4f}i", color=INK,
                 fontsize=12, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def plot_spacing(gamma: np.ndarray, bins: int = 60) -> Figure:
    import matplotlib.pyplot as plt

    gaps = normalized_spacings(gamma)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURFACE)
    s = np.linspace(0, 3.5, 400)
    left.hist(gaps, bins=bins, range=(0, 3.5), density=True, color=ORDINAL[0], edgecolor=SURFACE,
              linewidth=0.5, label=f"{len(gaps):,} zero spacings")
    left.plot(s, gue_spacing(s), color=INK, linewidth=1.8, label="GUE (sine kernel)")
    left.plot(s, wigner_surmise(s), color=ACCENT, linewidth=1.2, linestyle="--",
              label="Wigner surmise")
    left.plot(s, np.exp(-s), color=INK_SECONDARY, linewidth=1.2, linestyle=":", label="Poisson")
    left.set_xlim(0, 3.5)
    left.set_xlabel("normalized spacing", color=INK)
    left.set_ylabel("density", color=INK)
    left.set_title("Nearest-neighbor spacings", color=INK, fontsize=11, loc="left")

    u, r = pair_correlation(unfold(gamma), u_max=3.0, bins=bins)
    uu = np.linspace(0, 3, 400)
    right.bar(u, r, width=u[1] - u[0], color=ORDINAL[0], edgecolor=SURFACE, linewidth=0.5,
              label="zeros")
    right.plot(uu, montgomery(uu), color=INK, linewidth=1.8, label="1 − (sin πu / πu)²")
    right.axhline(1, color=INK_SECONDARY, linewidth=1.2, linestyle=":", label="Poisson")
    right.set_xlim(0, 3)
    right.set_xlabel("u (normalized distance)", color=INK)
    right.set_title("Pair correlation (Montgomery)", color=INK, fontsize=11, loc="left")

    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    left.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper right")
    right.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="lower right")
    fig.suptitle(
        f"Zeros {gamma[0]:,.0f} < t < {gamma[-1]:,.0f} against random matrix theory",
        color=INK, fontsize=12, x=0.01, ha="left",
    )
    fig.tight_layout()
    return fig


def plot_explicit(
    gamma: np.ndarray, counts: tuple[int, ...] = (10, 100, 1000), x_max: float = 100.0,
    samples: int = 3000,
) -> Figure:
    import matplotlib.pyplot as plt

    x = np.linspace(2.0, x_max, samples)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.8), facecolor=SURFACE)
    left.step(x, prime_pi(x), where="post", color=INK, linewidth=1.4, label="π(x)")
    left.plot(x, riemann_r(x), color=INK_SECONDARY, linewidth=1.2, linestyle="--",
              label="R(x), no zeros")
    right.step(x, chebyshev_psi(x) - x, where="post", color=INK, linewidth=1.4, label="ψ(x) − x")
    for K, color in zip(counts, (ORDINAL[0], ORDINAL[2], ACCENT), strict=True):
        label = f"{K:,} zeros"
        left.plot(x, pi_from_zeros(x, gamma[:K]), color=color, linewidth=1.2, label=label)
        right.plot(x, psi_from_zeros(x, gamma[:K]) - x, color=color, linewidth=1.2, label=label)
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
        ax.set_xlim(2, x_max)
        ax.set_xlabel("x", color=INK)
    left.set_title("Prime counting function", color=INK, fontsize=11, loc="left")
    right.set_title("Chebyshev ψ(x) minus its main term", color=INK, fontsize=11, loc="left")
    left.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper left")
    right.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="lower left", ncol=2)
    fig.suptitle("Rebuilding the primes from the zeros of ζ (explicit formulas, truncated)",
                 color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def plot_lehmer(n: int = 6709, half_width: float = 10.0, samples: int = 2000) -> Figure:
    import matplotlib.pyplot as plt

    gamma = first_zeros(n + 64)
    a, b = gamma[n - 1], gamma[n]
    t_min, t_max = a - half_width, a + half_width
    pad = 2 * (b - a)
    t = np.linspace(t_min, t_max, samples)
    zoom = np.linspace(a - pad, b + pad, samples // 2)

    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURFACE,
                                      gridspec_kw={"width_ratios": [1.6, 1]})
    origin = round(a)
    for ax, ts, shift in ((left, t, 0), (right, zoom, origin)):
        ax.set_facecolor(SURFACE)
        ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
        ax.plot(ts - shift, z_values(ts), color=SERIES, linewidth=2)
        _style(ax)
    left.set_xlabel("t", color=INK)
    right.set_xlabel(f"t − {origin}", color=INK)
    inside = gamma[(gamma >= t_min) & (gamma <= t_max)]
    left.plot(inside, np.zeros_like(inside), "o", color=SERIES, markersize=6,
              markeredgecolor=SURFACE, markeredgewidth=1.5)
    left.axvspan(a - pad, b + pad, color=ACCENT, alpha=0.15, linewidth=0)
    left.set_xlim(t_min, t_max)
    left.set_ylabel("Z(t)", color=INK)
    left.set_title(f"Z(t) near t = {a:.0f}", color=INK, fontsize=11, loc="left")
    right.plot([a - origin, b - origin], [0, 0], "o", color=ACCENT, markersize=7,
               markeredgecolor=SURFACE, markeredgewidth=1.5)
    right.set_xlim(zoom[0] - origin, zoom[-1] - origin)
    right.set_title(f"Zeros {n} and {n + 1}, {b - a:.4f} apart", color=INK, fontsize=11,
                    loc="left")
    fig.suptitle("A Lehmer pair: Z barely crosses zero twice", color=INK, fontsize=12, x=0.01,
                 ha="left")
    fig.tight_layout()
    return fig


def plot_log_z(sample: np.ndarray, T: float, N: int, bins: int = 120) -> Figure:
    import matplotlib.pyplot as plt

    x = np.linspace(-8, 5, 600)
    k = sample_cumulants(sample)
    selberg = 0.5 * np.log(np.log(T))
    cue = cue_density(x, N)
    hybrid = hybrid_density(x, N)
    hybrid_var = cue_cumulants(N)[1] + arithmetic_cumulants()[1]
    fig, (ax, diff) = plt.subplots(1, 2, figsize=(13, 4.8), facecolor=SURFACE,
                                   gridspec_kw={"width_ratios": [1.5, 1]})
    for a in (ax, diff):
        a.set_facecolor(SURFACE)
        _style(a)
        a.set_xlim(-8, 5)
        a.set_xlabel("log |ζ(1/2 + it)|", color=INK)
    label = f"log|ζ(1/2 + it)|, {len(sample):,} random t, variance {k[1]:.3f}"
    hist, edges, _ = ax.hist(sample, bins=bins, range=(-8, 5), density=True, color=ORDINAL[0],
                             edgecolor=SURFACE, linewidth=0.4, label=label)
    ax.plot(x, cue, color=INK, linewidth=1.8,
            label=f"CUE, N = {N}, variance {cue_cumulants(N)[1]:.3f}")
    ax.plot(x, hybrid, color=ORDINAL[4], linewidth=1.4, linestyle=":",
            label=f"CUE with the arithmetic factor, variance {hybrid_var:.3f}")
    ax.plot(x, np.exp(-(x**2) / (2 * selberg)) / np.sqrt(2 * np.pi * selberg), color=ACCENT,
            linewidth=1.4, linestyle="--", label=f"Selberg normal, variance {selberg:.3f}")
    ax.set_ylabel("density", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper left")
    ax.set_title(f"Value distribution near t = {T:.3g} against random unitary matrices",
                 color=INK, fontsize=12, loc="left")
    mid = (edges[:-1] + edges[1:]) / 2
    err = np.sqrt(hist / (len(sample) * (edges[1] - edges[0])))
    diff.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    diff.fill_between(mid, -2 * err, 2 * err, color=GRID, step="mid", label="±2σ sampling noise")
    diff.plot(mid, hist - cue_density(mid, N), color=ORDINAL[2], linewidth=1.2,
              drawstyle="steps-mid", label="histogram − CUE")
    diff.plot(x, hybrid - cue, color=ACCENT, linewidth=1.8, label="arithmetic factor − CUE")
    diff.set_ylabel("density difference", color=INK)
    diff.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="lower left")
    diff.set_title("What CUE misses", color=INK, fontsize=12, loc="left")
    fig.tight_layout()
    return fig


def plot_l_z(
    labels: tuple[tuple[int, int], ...] = ((4, 3), (5, 2), (7, 3)), t_min: float = -30.0,
    t_max: float = 30.0, samples: int = 2000,
) -> Figure:
    import matplotlib.pyplot as plt

    t = np.linspace(t_min, t_max, samples)
    fig, axes = plt.subplots(len(labels), 1, figsize=(10, 2.3 * len(labels) + 0.6), sharex=True,
                             facecolor=SURFACE)
    for ax, (q, n) in zip(np.atleast_1d(axes), labels, strict=True):
        chi = conrey(q, n)
        zeros = l_zeros(chi, t_min, t_max)
        ax.set_facecolor(SURFACE)
        ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
        ax.plot(t, l_z(t, chi), color=SERIES, linewidth=1.6)
        ax.plot(zeros, np.zeros_like(zeros), "o", color=SERIES, markersize=6,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
        _style(ax)
        kind = "real" if chi.is_real else f"order {chi.order}"
        parity = "even" if chi.parity == 0 else "odd"
        ax.set_title(f"χ_{q}({n}, ·): {kind}, {parity}, ε = {_fmt_unit(chi.root_number())}",
                     color=INK, fontsize=10, loc="left")
        ax.set_ylabel(f"Z(t, χ_{q}({n},·))", color=INK, fontsize=9)
    ax.set_xlim(t_min, t_max)
    ax.set_xlabel("t", color=INK)
    fig.suptitle("Hardy Z functions of Dirichlet L-functions: zeros of L(1/2 + it, χ)",
                 color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def _fmt_unit(z: complex) -> str:
    if abs(z.imag) < 1e-12:
        return f"{z.real:+.0f}"
    return f"e^({np.angle(z) / (2 * np.pi):.4f}·2πi)"


def plot_l_zeros(q: int = 31, T: float = 20.0) -> Figure:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(11, 5.4), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    chars = dict(primitive_characters(q))
    order = []
    for n in chars:
        if n not in order:
            order += [n] if pow(n, -1, q) == n else [n, pow(n, -1, q)]
    seen = set()
    for i, n in enumerate(order):
        chi = chars[n]
        zeros = l_zeros(chi, -T, T)
        if chi.is_real:
            color, label = ACCENT, "real character"
        elif chi.parity == 0:
            color, label = ORDINAL[1], "even, complex"
        else:
            color, label = ORDINAL[4], "odd, complex"
        ax.plot(np.full_like(zeros, i), zeros, "_", color=color, markersize=9, markeredgewidth=2,
                label=None if label in seen else label)
        seen.add(label)
    ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    _style(ax)
    ax.set_xlim(-1, len(order))
    ax.set_ylim(-T, T)
    ax.set_xticks(range(len(order)), [str(n) for n in order], fontsize=8)
    ax.set_xlabel(f"Conrey label n of χ_{q}(n, ·), next to its conjugate n⁻¹ mod {q}", color=INK)
    ax.set_ylabel("γ, where L(1/2 + iγ, χ) = 0", color=INK)
    ax.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper left",
              bbox_to_anchor=(1.0, 1.0))
    ax.set_title(f"Zeros of the primitive L-functions mod {q}: χ and its conjugate mirror each "
                 "other, real characters are symmetric", color=INK, fontsize=11, loc="left")
    fig.tight_layout()
    return fig


def plot_low_lying(
    symplectic: OneLevel, unitary: OneLevel, labels: tuple[str, str], bins: int = 30,
    x_max: float = 3.0,
) -> Figure:
    import matplotlib.pyplot as plt

    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.8), facecolor=SURFACE)
    x = np.linspace(0, x_max, 400)
    for fam, label, color in ((unitary, labels[1], ORDINAL[0]), (symplectic, labels[0], SERIES)):
        c, d = fam.density(bins, x_max)
        left.step(c, d, where="mid", color=color, linewidth=1.8, label=label)
    left.plot(x, w_symplectic(x), color=INK, linewidth=1.4, label="USp: 1 − sin 2πx / 2πx")
    left.plot(x, w_unitary(x), color=INK_SECONDARY, linewidth=1.2, linestyle=":", label="U: 1")
    left.set_xlim(0, x_max)
    left.set_ylim(0, 2.0)
    left.set_xlabel("x (zero height in units of the mean spacing)", color=INK)
    left.set_ylabel("zeros per L-function per unit x", color=INK)
    left.set_title("One-level density near s = 1/2", color=INK, fontsize=11, loc="left")

    for fam, label, color in ((unitary, labels[1], ORDINAL[0]), (symplectic, labels[0], SERIES)):
        right.hist(fam.first(), bins=bins, range=(0, x_max / 2), density=True, histtype="step",
                   color=color, linewidth=1.8, label=label)
    right.set_xlim(0, x_max / 2)
    right.set_ylim(0, 2.0)
    right.set_xlabel("x of the lowest zero", color=INK)
    right.set_ylabel("density", color=INK)
    right.set_title("Where the first zero sits", color=INK, fontsize=11, loc="left")

    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    left.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper center", ncol=2)
    right.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper right")
    fig.suptitle("Low-lying zeros of Dirichlet L-functions against Katz-Sarnak",
                 color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout()
    return fig



def plot_low_lying_scales(
    families: dict[float, OneLevel], window: tuple[float, float] = (0.0, 0.3), bins: int = 30,
    x_max: float = 3.0, ratios: dict[float, np.ndarray] | None = None,
) -> Figure:
    import matplotlib.pyplot as plt

    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURFACE,
                                      gridspec_kw={"width_ratios": [1.6, 1]})
    x = np.linspace(0, x_max, 400)
    colors = ORDINAL[::2] if len(families) <= 3 else ORDINAL
    lo, hi = window
    window_x = np.linspace(lo, hi, 2001)
    expected = np.mean(w_symplectic(window_x))
    size, excess, err = [], [], []
    for (d, fam), color in zip(sorted(families.items()), colors[: len(families)], strict=True):
        c, dens = fam.density(bins, x_max)
        left.step(c, dens, where="mid", color=color, linewidth=1.8,
                  label=f"|d| ≈ {_fmt_power(d)}, {len(fam.zeros):,} characters")
        k = np.count_nonzero((fam.x > lo) & (fam.x <= hi))
        scale = len(fam.zeros) * (hi - lo)
        size.append(d)
        excess.append(k / scale - expected)
        err.append(np.sqrt(k) / scale)
        right.errorbar(1 / np.log(d), excess[-1], yerr=err[-1], fmt="o", color=color,
                       markersize=7, capsize=3, linewidth=1.6)
        if ratios is not None:
            rho = quadratic_density(ratios[d], np.concatenate([x, window_x]))
            left.plot(x, rho[: len(x)], color=color, linewidth=1.2, linestyle="--")
            right.plot(1 / np.log(d), np.mean(rho[len(x) :]) - expected, "D", color=color,
                       markerfacecolor=SURFACE, markersize=7, markeredgewidth=1.6)
    left.plot(x, w_symplectic(x), color=INK, linewidth=1.4, label="USp: 1 − sin 2πx / 2πx")
    if ratios is not None:
        left.plot([], [], color=INK_SECONDARY, linewidth=1.2, linestyle="--",
                  label="ratios conjecture, same characters")
    left.set_xlim(0, x_max)
    left.set_ylim(0, 2.2)
    left.set_xlabel("x (zero height in units of the mean spacing)", color=INK)
    left.set_ylabel("zeros per L-function per unit x", color=INK)
    left.set_title(f"Quadratic characters at {len(families)} conductor sizes", color=INK,
                   fontsize=11, loc="left")
    left.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper center", ncol=2)

    u = 1 / np.log(size)
    wts = 1 / np.square(err)
    slope = np.sum(wts * u * np.asarray(excess)) / np.sum(wts * u * u)
    grid = np.linspace(0, u.max() * 1.15, 50)
    right.plot(grid, slope * grid, color=INK_SECONDARY, linewidth=1.2, linestyle="--",
               label=f"{slope:.2f} / log |d|")
    if ratios is not None:
        big = np.logspace(np.log10(min(size)), 60, 14)
        curve = [np.mean(quadratic_density([-D, D], window_x, n=1000)) - expected for D in big]
        right.plot(1 / np.log(big), curve, color=ACCENT, linewidth=1.4,
                   label="ratios conjecture, d = ±D")
        right.plot([], [], "D", color=INK_SECONDARY, markerfacecolor=SURFACE, markersize=7,
                   markeredgewidth=1.6, label="ratios conjecture, same characters")
    right.axhline(0, color=INK_SECONDARY, linewidth=1, linestyle=":")
    right.set_xlim(0, grid[-1])
    right.set_ylim(min(0, min(e - s for e, s in zip(excess, err, strict=True)) * 1.2),
                   max(e + s for e, s in zip(excess, err, strict=True)) * 1.25)
    right.set_xlabel("1 / log |d|", color=INK)
    right.set_ylabel(f"density − W_USp on ({lo:g}, {hi:g}]", color=INK)
    right.set_title("Excess near the central point", color=INK, fontsize=11, loc="left")
    right.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper left")
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_first_zero(
    families: dict[float, OneLevel], n_grid: ArrayLike | None = None, s_max: float = 1.5,
    bins: int = 30,
) -> Figure:
    import matplotlib.pyplot as plt

    grid = np.arange(1.4, 6.01, 0.05) if n_grid is None else np.asarray(n_grid, dtype=float)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURFACE)
    s = np.linspace(0, s_max, 301)
    colors = ORDINAL[::2] if len(families) <= 3 else ORDINAL
    for (d, fam), color in zip(sorted(families.items()), colors[: len(families)], strict=True):
        first = np.sort(fam.first())
        n, lo, hi, _ = fit_n(first, grid)
        left.hist(first, bins=bins, range=(0, s_max), density=True, histtype="step",
                  color=color, linewidth=1.8, label=f"|d| ≈ {_fmt_power(d)}")
        left.plot(s, first_density(symplectic_n(n), s), color=color, linewidth=1.2,
                  linestyle="--", label=f"USp(2N), N = {n:.2f}")
        survival = 1 - np.arange(len(first)) / len(first)
        right.step(first, survival, where="post", color=color, linewidth=1.8,
                   label=f"|d| ≈ {_fmt_power(d)}, N = {n:.2f} ({lo:.2f} to {hi:.2f})")
        right.plot(s, gap(symplectic_n(n), s), color=color, linewidth=1.2, linestyle="--")
    left.plot(s, first_density(symplectic, s), color=INK, linewidth=1.4, label="USp, N → ∞")
    left.plot(s, first_density(even_orthogonal, s), color=INK_SECONDARY, linewidth=1.2,
              linestyle=":", label="SO(even), N → ∞")
    right.plot(s, gap(symplectic, s), color=INK, linewidth=1.4, label="USp, N → ∞")
    right.plot(s, gap(even_orthogonal, s), color=INK_SECONDARY, linewidth=1.2, linestyle=":",
               label="SO(even), N → ∞")
    left.set_xlim(0, s_max)
    left.set_ylim(0, 2.6)
    left.set_xlabel("x of the lowest zero", color=INK)
    left.set_ylabel("density", color=INK)
    left.set_title("Lowest zero of L(s, χ_d) against Fredholm determinants", color=INK,
                   fontsize=11, loc="left")
    left.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper right", ncol=2)
    right.set_yscale("log")
    right.set_xlim(0, s_max)
    right.set_ylim(1e-3, 1.2)
    right.set_xlabel("s", color=INK)
    right.set_ylabel("P(lowest zero > s)", color=INK)
    right.set_title("Gap probability, dashed: fitted USp(2N)", color=INK, fontsize=11,
                    loc="left")
    right.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="lower left")
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig

def _fmt_power(x: float) -> str:
    e = int(np.floor(np.log10(x)))
    m = x / 10**e
    sup = str(e).translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"10{sup}" if round(m, 1) == 1 else f"{m:.1f}·10{sup}"


def plot_landscape(
    t_min: float = 0.0, t_max: float = 40.0, sigma_min: float = -0.5, sigma_max: float = 2.0,
    floor: float = -1.5, ceiling: float = 1.0, width: int = 600,
) -> Figure:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LightSource

    zeros = find_zeros(max(t_min, 10.0), t_max)
    t = np.union1d(np.linspace(t_min, t_max, width), zeros)
    rows = int(round(width * (sigma_max - sigma_min) / (t_max - t_min))) * 2
    sigma = np.union1d(np.linspace(sigma_min, sigma_max, rows), [0.5, 1.0])
    S, T = np.meshgrid(sigma, t, indexing="ij")
    s = S + 1j * T
    w = np.full(s.shape, ceiling)
    finite = s != 1
    with np.errstate(divide="ignore"):
        w[finite] = np.clip(np.log10(np.abs(zeta(s[finite]))), floor, ceiling)

    fig = plt.figure(figsize=(11, 6), facecolor=SURFACE)
    ax = fig.add_subplot(projection="3d", computed_zorder=False)
    ax.set_facecolor(SURFACE)
    light = LightSource(azdeg=300, altdeg=40)
    rgb = light.shade(w, cmap=plt.get_cmap("viridis"), vmin=floor, vmax=ceiling,
                      blend_mode="soft", vert_exag=0.3)
    ax.plot_surface(T, S, w, facecolors=rgb, rstride=1, cstride=1, linewidth=0,
                    antialiased=False, shade=False)
    ax.scatter(zeros, np.full(len(zeros), 0.5), np.full(len(zeros), floor), color=ACCENT, s=16,
               depthshade=False)
    ax.view_init(elev=42, azim=-70)
    ax.set_box_aspect((2.6, 1.1, 0.7), zoom=1.2)
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(sigma_min, sigma_max)
    ax.set_zlim(floor, ceiling)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color(SURFACE)
        axis._axinfo["grid"]["color"] = GRID
    ax.tick_params(colors=INK_SECONDARY, labelsize=8)
    ax.set_xlabel("t = Im s", color=INK)
    ax.set_ylabel("σ = Re s", color=INK)
    ax.set_zlabel("log₁₀ |ζ(s)|", color=INK)
    fig.text(0.02, 0.95, f"log |ζ(σ + it)|: a spike at the pole s = 1 and a pit at each of the "
             f"{len(zeros)} zeros on σ = 1/2", color=INK, fontsize=12)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    return fig


def sweep_sigma(
    t_max: float = 40.0, frames: int = 60, samples: int = 1500, amplitude: float = 0.4,
):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    t = np.linspace(0, t_max, samples)
    sigmas = 0.5 + amplitude * np.sin(2 * np.pi * np.arange(frames) / frames) ** 3
    curves = zeta(sigmas[:, None] + 1j * t[None, :])
    zeros = find_zeros(10, t_max)

    fig, (left, right) = plt.subplots(1, 2, figsize=(10, 4.6), facecolor=SURFACE,
                                      gridspec_kw={"width_ratios": [1, 1.35]})
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    reach = np.abs(curves[:, t > 5]).max() * 1.05
    spiral, = left.plot([], [], color=SERIES, linewidth=1.2)
    left.plot([0], [0], "o", color=ACCENT, markersize=6, markeredgecolor=SURFACE)
    left.set_xlim(-reach, reach)
    left.set_ylim(-reach, reach)
    left.set_aspect("equal")
    left.set_xlabel("Re ζ(σ + it)", color=INK)
    left.set_ylabel("Im ζ(σ + it)", color=INK)
    modulus, = right.plot([], [], color=SERIES, linewidth=1.4)
    right.plot(zeros, np.zeros(len(zeros)), "o", color=ACCENT, markersize=6,
               markeredgecolor=SURFACE)
    right.set_xlim(0, t_max)
    right.set_ylim(0, reach)
    right.set_xlabel("t", color=INK)
    right.set_ylabel("|ζ(σ + it)|", color=INK)
    title = fig.suptitle("", color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout()

    def draw(i):
        w = curves[i]
        spiral.set_data(w.real, w.imag)
        modulus.set_data(t, np.abs(w))
        title.set_text(f"σ = {sigmas[i]:.2f}: ζ(σ + it) for 0 ≤ t ≤ {t_max:g}, "
                       "zeros only when σ = 1/2")
        return spiral, modulus, title

    return FuncAnimation(fig, draw, frames=frames, interval=80, blit=False)


def plot_multieval(heights: ArrayLike = (1e4, 1e6, 1e8, 1e10, 1e12), points: int = 4) -> Figure:
    import time

    import matplotlib.pyplot as plt
    import mpmath

    heights = np.asarray(heights, dtype=float)
    err_direct, err_window, us_direct, us_window = [], [], [], []
    for T in heights:
        win = ZWindow(T, T + 100)
        u = np.linspace(win.u_lo, win.u_hi, points + 2)[1:-1]
        with mpmath.workdps(35):
            exact = np.array([float(mpmath.siegelz(mpmath.mpf(win.t0) + mpmath.mpf(x)))
                              for x in u])
        err_direct.append(np.abs(z_rs(win.t0 + u) - exact).max())
        err_window.append(np.abs(win.z(u) - exact).max())
        dense = np.linspace(win.u_lo, win.u_hi, 20_000)
        start = time.perf_counter()
        win.z(dense)
        us_window.append((time.perf_counter() - start) / len(dense) * 1e6)
        few = dense[:200]
        start = time.perf_counter()
        z_rs(win.t0 + few)
        us_direct.append((time.perf_counter() - start) / len(few) * 1e6)

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.4), facecolor=SURFACE)
    for ax, direct, window, label in (
        (left, err_direct, err_window, "max |Z − mpmath|"),
        (right, us_direct, us_window, "μs per point"),
    ):
        ax.set_facecolor(SURFACE)
        ax.loglog(heights, direct, "o-", color=ORDINAL[1], linewidth=2, label="direct sum")
        ax.loglog(heights, window, "o-", color=ACCENT, linewidth=2, label="NUFFT window")
        _style(ax)
        ax.set_xlabel("t", color=INK)
        ax.set_ylabel(label, color=INK)
    left.set_title("Accuracy against mpmath", color=INK, fontsize=11, loc="left")
    right.set_title("Cost per evaluation after the window is built", color=INK, fontsize=11,
                    loc="left")
    left.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper left")
    fig.tight_layout()
    return fig



def plot_afe(
    conductors: ArrayLike = (101, 1009, 10007, 100003), heights: ArrayLike = (1, 10, 100, 1000),
    q_fixed: int = 10007, t_fixed: float = 10.0, points: int = 8,
) -> Figure:
    import time

    import matplotlib.pyplot as plt

    def ms(f, s):
        start = time.perf_counter()
        f(s)
        return (time.perf_counter() - start) / len(s) * 1e3

    fig, (left, right) = plt.subplots(1, 2, figsize=(11, 4.4), facecolor=SURFACE)
    for ax, x, pairs, label in (
        (left, conductors, [(q, t_fixed) for q in conductors], "conductor q"),
        (right, heights, [(q_fixed, t) for t in heights], "t"),
    ):
        em, smooth = [], []
        for q, t in pairs:
            chi = conrey(q, 2)
            chi.root_number()
            s = 0.5 + 1j * np.linspace(t, t + 1, points)
            em.append(ms(lambda s, chi=chi: dirichlet_l(s, chi, method="em"), s))
            smooth.append(ms(lambda s, chi=chi: l_afe(s, chi), s))
        ax.set_facecolor(SURFACE)
        ax.loglog(x, em, "o-", color=ORDINAL[1], linewidth=2, label="Euler-Maclaurin")
        ax.loglog(x, smooth, "o-", color=ACCENT, linewidth=2, label="smoothed AFE")
        for xi, yi, (q, t) in zip(x, smooth, pairs, strict=True):
            ax.annotate(f"{2 * afe_length(q, t)} terms", (xi, yi), textcoords="offset points",
                        xytext=(0, 8), ha="center", fontsize=8, color=INK_SECONDARY)
        _style(ax)
        ax.set_xlabel(label, color=INK)
        ax.set_ylabel("ms per value of L(1/2 + it, χ)", color=INK)
    left.set_title(f"Cost against the conductor at t = {t_fixed:g}", color=INK, fontsize=11,
                   loc="left")
    right.set_title(f"Cost against the height for q = {q_fixed}", color=INK, fontsize=11,
                    loc="left")
    left.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper left")
    fig.tight_layout()
    return fig


CURVES = {
    "11a3": Curve(0, -1, 1, 0, 0),
    "37a1": Curve(0, 0, 1, -1, 0),
    "389a1": Curve(0, 1, 1, -2, 0),
    "5077a1": Curve(0, 0, 1, -7, 6),
}


def plot_ec_z(t_max: float = 6.0, samples: int = 601) -> Figure:
    import matplotlib.pyplot as plt

    t = np.linspace(-t_max, t_max, samples)
    fig, axes = plt.subplots(len(CURVES), 1, figsize=(10, 2.1 * len(CURVES) + 0.6), sharex=True,
                             facecolor=SURFACE)
    for ax, (name, E) in zip(axes, CURVES.items(), strict=True):
        L = EllipticL(E)
        z = zeros_below(L, t_max)
        ax.set_facecolor(SURFACE)
        ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
        ax.plot(t, hardy_z(t, L), color=SERIES, linewidth=1.6)
        zeros = np.concatenate([-z.t, z.t])
        ax.plot(zeros, np.zeros_like(zeros), "o", color=SERIES, markersize=5,
                markeredgecolor=SURFACE, markeredgewidth=1.2)
        ax.plot([0], [0], "o", color=ACCENT, markersize=8 if z.rank else 0,
                markeredgecolor=SURFACE, markeredgewidth=1.2)
        _style(ax)
        w = "+1" if L.root_number == 1 else "−1"
        ax.set_title(f"{name}: {_fmt_curve(E)},  N = {L.conductor},  w = {w},  "
                     f"order {z.rank} at s = 1/2", color=INK, fontsize=10, loc="left")
        ax.set_ylabel("Z(t)", color=INK, fontsize=9)
    ax.set_xlim(-t_max, t_max)
    ax.set_xlabel("t", color=INK)
    fig.suptitle("Hardy Z of elliptic curve L-functions: the central zero has order 0, 1, 2, 3",
                 color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def _fmt_curve(E: Curve) -> str:
    lhs = "y²" + "".join(f" + {m}" for m, c in (("xy", E.a1), ("y", E.a3)) if c)
    rhs = "".join(f" {'+' if c > 0 else '−'} {abs(c) if abs(c) != 1 or not m else ''}{m}"
                  for m, c in (("x²", E.a2), ("x", E.a4), ("", E.a6)) if c)
    return f"{lhs} = x³{rhs}"


def plot_central_values(d: ArrayLike, value: ArrayLike, kappa: float, c_max: int = 6) -> Figure:
    import matplotlib.pyplot as plt

    d = np.asarray(d, dtype=float)
    v = np.asarray(value, dtype=float)
    fig, ax = plt.subplots(figsize=(11, 5), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)
    grid = np.linspace(d.min(), d.max(), 400)
    for c in range(1, c_max + 1):
        ax.plot(grid, kappa * c * c / np.sqrt(grid), color=GRID, linewidth=1.2, zorder=1)
        ax.annotate(f"c = {c}", (grid[-1], kappa * c * c / np.sqrt(grid[-1])),
                    xytext=(4, 0), textcoords="offset points", va="center", fontsize=8,
                    color=INK_SECONDARY)
    zero = np.abs(v) < 1e-6
    ax.plot(d[~zero], v[~zero], ".", color=SERIES, markersize=4, zorder=2,
            label=f"L(E_d, 1/2) ≠ 0 ({np.count_nonzero(~zero)})")
    ax.plot(d[zero], v[zero], "|", color=ACCENT, markersize=7, zorder=2,
            label=f"L(E_d, 1/2) = 0, rank 2 ({np.count_nonzero(zero)})")
    _style(ax)
    ax.set_xlim(0, d.max() * 1.07)
    ax.set_ylim(-0.08, kappa * c_max * c_max / np.sqrt(d.max()) * 1.4)
    ax.set_xlabel("fundamental discriminant d > 0 with w(E_d) = +1", color=INK)
    ax.set_ylabel("L(E_d, 1/2)", color=INK)
    ax.set_title(f"Central values of the even twists of 11a3 sit on κ c² / √d, κ = {kappa}, "
                 "c an integer (Waldspurger, Kohnen-Zagier)", color=INK, fontsize=11, loc="left")
    ax.legend(frameon=True, facecolor=SURFACE, edgecolor=SURFACE, framealpha=1, fontsize=9,
              labelcolor=INK_SECONDARY, loc="upper right")
    fig.tight_layout()
    return fig

def plot_twist_zeros(
    even: ArrayLike, odd: ArrayLike, excised: ArrayLike, n: int, cutoff: float,
    n_grid: ArrayLike | None = None, s_max: float = 2.0, bins: int = 32,
) -> Figure:
    import matplotlib.pyplot as plt

    grid = np.arange(1.0, 8.01, 0.05) if n_grid is None else np.asarray(n_grid, dtype=float)
    even, odd = np.sort(even), np.sort(odd)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.8), facecolor=SURFACE)
    s = np.linspace(0, s_max, 301)
    for ax, x in ((left, even), (right, odd)):
        ax.hist(x, bins=bins, range=(0, s_max), density=True, histtype="stepfilled",
                color=ORDINAL[0], alpha=0.55, label=f"twists ({len(x)})")
    n_fit, lo, hi, _ = fit_n(even, grid, even_orthogonal_n)
    left.plot(s, first_density(even_orthogonal, s), color=INK_SECONDARY, linewidth=1.3,
              linestyle=":", label="SO(even), N → ∞")
    left.plot(s, first_density(even_orthogonal_n(n_fit), s), color=SERIES, linewidth=1.6,
              linestyle="--", label=f"SO(2N), fitted N = {n_fit:.2f} ({lo:.2f} to {hi:.2f})")
    left.hist(np.asarray(excised, dtype=float), bins=bins, range=(0, s_max), density=True,
              histtype="step", color=ACCENT, linewidth=1.8,
              label=f"SO({2 * n}) with Λ(1) ≥ {cutoff:.4f}")
    right.plot(s, first_density(symplectic, s), color=INK_SECONDARY, linewidth=1.3,
               linestyle=":", label="SO(odd), N → ∞")
    right.plot(s, first_density(odd_orthogonal_n(n), s), color=SERIES, linewidth=1.6,
               linestyle="--", label=f"SO({2 * n + 1})")
    left.set_ylabel("density", color=INK)
    left.set_title("Even twists, w = +1, rank 0", color=INK, fontsize=11, loc="left")
    right.set_title("Odd twists, w = −1, rank 1: first zero after the central one",
                    color=INK, fontsize=11, loc="left")
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        ax.set_xlim(0, s_max)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.3)
        ax.set_xlabel("x of the lowest zero above s = 1/2", color=INK)
        ax.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
        _style(ax)
    fig.tight_layout()
    return fig


def plot_rank_two(
    first: ArrayLike, size: ArrayLike, n: int, sample: ArrayLike, sample_size: ArrayLike,
    s_max: float = 3.0, bins: int = 24, excised: ArrayLike | None = None, removed: float = 0.0,
) -> Figure:
    import matplotlib.pyplot as plt

    x = np.asarray(first, dtype=float)
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(12, 4.8), facecolor=SURFACE, gridspec_kw={"width_ratios": [1.3, 1]}
    )
    s = np.linspace(0, s_max, 241)
    left.hist(x, bins=bins, range=(0, s_max), density=True, histtype="stepfilled",
              color=ORDINAL[0], alpha=0.55, label=f"rank 2 twists ({len(x)})")
    left.plot(s, first_density(even_orthogonal, s), color=INK_SECONDARY, linewidth=1.2,
              linestyle=":", label="SO(even), no forced pair")
    curves = (
        (independent_n(n), ORDINAL[2], "--", f"independent: SO({2 * n - 2})"),
        (conditioned_n(n), ORDINAL[4], "-.", "Λ(1) → 0 inside SO(2N)"),
        (interaction_n(n), ACCENT, "-", "interaction: pair at 1"),
    )
    for k, color, style, label in curves:
        left.plot(s, first_density(k, s), color=color, linewidth=1.7, linestyle=style,
                  label=f"{label}, N = {n}")
    left.plot(s, first_density(bessel(1.5), s), color=ACCENT, linewidth=1.1, linestyle=":",
              label="interaction, N → ∞")
    if excised is not None:
        left.hist(np.asarray(excised, dtype=float), bins=2 * bins, range=(0, s_max), density=True,
                  histtype="step", color=INK, linewidth=1.3,
                  label=f"interaction, smallest {removed:.0%} of Λ''(1) removed")
    left.set_xlim(0, s_max)
    left.set_ylim(0, left.get_ylim()[1] * 1.25)
    left.set_xlabel("x of the lowest zero above the double zero at s = 1/2", color=INK)
    left.set_ylabel("density", color=INK)
    left.set_title("Rank 2 even twists of 11a3 against SO(2N) with a pair at 1", color=INK,
                   fontsize=11, loc="left")
    left.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
    right.scatter(np.log(np.asarray(sample_size, dtype=float)), np.asarray(sample, dtype=float),
                  s=4, color=GRID, alpha=0.6, linewidths=0, label=f"interaction model, N = {n}")
    right.scatter(np.log(np.asarray(size, dtype=float)), x, s=10, color=SERIES, alpha=0.8,
                  linewidths=0, label="rank 2 twists")
    right.set_ylim(0, s_max)
    right.set_xlabel("log of the leading Taylor coefficient (normalized)", color=INK)
    right.set_ylabel("x of the lowest zero", color=INK)
    right.set_title("Small second derivative, close first zero", color=INK, fontsize=11,
                    loc="left")
    right.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper left",
                 markerscale=3)
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_odd_excision(
    first: ArrayLike, theta: np.ndarray, predicted: dict[str, float], c_max: float = 300.0,
    s_max: float = 1.8, bins: int = 18,
) -> Figure:
    import matplotlib.pyplot as plt

    x = np.sort(np.asarray(first, dtype=float))
    n = theta.shape[1]
    lam = char_poly_at_one(theta) * np.exp(n / 2)
    eig = n * theta[:, 0] / np.pi
    edges = np.linspace(0, s_max, bins + 1)
    cs = np.geomspace(0.1, c_max, 90)
    chi = np.array([binned_chi2(x, eig[lam >= c], edges) for c in cs])
    best = cs[np.argmin(chi)]
    kept = np.mean(lam >= best)
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(12, 4.8), facecolor=SURFACE, gridspec_kw={"width_ratios": [1.3, 1]}
    )
    left.hist(x, bins=edges, density=True, histtype="stepfilled", color=ORDINAL[0], alpha=0.55,
              label=f"odd twists, rank 1 ({len(x)})")
    left.hist(eig, bins=edges, density=True, histtype="step", color=INK_SECONDARY,
              linewidth=1.3, linestyle=":", label=f"SO({2 * n + 1})")
    name, c = max(predicted.items(), key=lambda kv: kv[1])
    left.hist(eig[lam >= c], bins=edges, density=True, histtype="step", color=SERIES,
              linewidth=1.6, label=f"SO({2 * n + 1}), |Λ'(1)| ≥ {c:.2f}·e^(−N/2) ({name})")
    left.hist(eig[lam >= best], bins=edges, density=True, histtype="step", color=ACCENT,
              linewidth=1.8, label=f"best cut c = {best:.0f}, keeps {kept:.0%}")
    left.set_xlim(0, s_max)
    left.set_ylim(0, left.get_ylim()[1] * 1.35)
    left.set_xlabel("x of the lowest zero above the central one", color=INK)
    left.set_ylabel("density", color=INK)
    left.set_title("Odd twists of 11a3 against excised SO(2N + 1)", color=INK, fontsize=11,
                   loc="left")
    right.plot(cs, chi, color=SERIES, linewidth=1.6)
    for (label, c), style in zip(predicted.items(), ("--", ":"), strict=False):
        right.axvline(c, color=INK_SECONDARY, linewidth=1.1, linestyle=style,
                      label=f"{label}: c = {c:.2f}")
    right.axvline(best, color=ACCENT, linewidth=1.1, label=f"minimum: c = {best:.0f}")
    right.set_xscale("log")
    right.set_xlim(cs[0], cs[-1])
    right.set_ylim(0, chi[0] * 1.3)
    right.set_xlabel("cut c in |Λ'(1)| ≥ c·e^(−N/2)", color=INK)
    right.set_ylabel(f"χ² on {bins} bins of (0, {s_max})", color=INK)
    right.set_title(f"N = {n}: the arithmetic cut does nothing", color=INK, fontsize=11,
                    loc="left")
    left.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
    right.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="lower left",
                 bbox_to_anchor=(0.42, 0.02))
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_twist_density(
    x: ArrayLike, families: dict[str, tuple[int, ArrayLike, int, ArrayLike]], x_max: float = 3.0,
    bins: int = 30,
) -> Figure:
    import matplotlib.pyplot as plt

    x = np.asarray(x, dtype=float)
    edges = np.linspace(0, x_max, bins + 1)
    mid = (edges[1:] + edges[:-1]) / 2
    width = edges[1] - edges[0]
    fig, axes = plt.subplots(1, len(families), figsize=(12, 4.8), facecolor=SURFACE)
    for ax, (title, (sign, zeros, count, pred)) in zip(axes, families.items(), strict=True):
        z = np.asarray(zeros, dtype=float)
        obs = np.histogram(z, edges)[0]
        limit = 1 + sign * np.sinc(2 * x)
        chi = {}
        for key, curve in (("ratios", pred), ("limit", limit)):
            expect = np.interp(mid, x, curve) * width * count
            chi[key] = float(np.sum((obs - expect) ** 2 / expect))
        ax.bar(mid, obs / (width * count), width=width, color=ORDINAL[0], alpha=0.55,
               label=f"{len(z)} zeros of {count} twists")
        ax.plot(x, pred, color=SERIES, linewidth=1.8,
                label=f"ratios conjecture (HKS), χ² = {chi['ratios']:.1f}")
        ax.plot(x, limit, color=INK_SECONDARY, linewidth=1.3, linestyle=":",
                label=f"{'SO(even)' if sign == 1 else 'SO(odd)'}, N → ∞, χ² = {chi['limit']:.1f}")
        ax.set_facecolor(SURFACE)
        ax.set_xlim(0, x_max)
        ax.set_ylim(0, max(2.2, ax.get_ylim()[1]) * 1.18)
        ax.set_xlabel("x = θ(t)/π, zeros above the center", color=INK)
        ax.set_title(title, color=INK, fontsize=11, loc="left")
        ax.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
        _style(ax)
    axes[0].set_ylabel(f"one-level density ({bins} bins, χ² over all of them)", color=INK)
    fig.tight_layout()
    return fig


def plot_odd_effective(
    first: ArrayLike, n_eff: ArrayLike, thetas: dict[int, np.ndarray], arithmetic: float,
    s_max: float = 1.8, bins: int = 18,
) -> Figure:
    import matplotlib.pyplot as plt

    x = np.sort(np.asarray(first, dtype=float))
    n_eff = np.asarray(n_eff, dtype=float)
    edges = np.linspace(0, s_max, bins + 1)
    drops = np.linspace(0, 0.12, 49)
    s = np.linspace(0, s_max, 181)
    mix = first_density_mix(n_eff, s, odd_orthogonal_half)
    expect = np.diff(np.interp(edges, s, np.r_[0, np.cumsum((mix[1:] + mix[:-1]) / 2
                                                         * np.diff(s))])) * len(x)
    chi_mix = float(np.sum((np.histogram(x, edges)[0] - expect) ** 2 / expect))
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(12, 4.8), facecolor=SURFACE, gridspec_kw={"width_ratios": [1.3, 1]}
    )
    left.hist(x, bins=edges, density=True, histtype="stepfilled", color=ORDINAL[0], alpha=0.55,
              label=f"odd twists, rank 1 ({len(x)})")
    left.plot(s, mix, color=SERIES, linewidth=1.8,
              label=f"SO(2N + 1), N_eff per twist (mean {n_eff.mean():.2f}), χ² = {chi_mix:.1f}")
    n0 = int(round(n_eff.mean()))
    styles = dict(zip(sorted(thetas), (":", "-", "--"), strict=False))
    for n, angles in sorted(thetas.items()):
        lam = char_poly_at_one(angles)
        eig = (n + 0.5) * angles[:, 0] / np.pi
        cuts = np.quantile(lam, drops)
        chi = np.array([binned_chi2(x, eig[lam >= c], edges) for c in cuts])
        right.plot(drops, chi, color=ACCENT if n == n0 else INK_SECONDARY,
                   linewidth=1.8 if n == n0 else 1.2, linestyle=styles[n],
                   label=f"SO({2 * n + 1}), best χ² = {chi.min():.1f}")
        if n == n0:
            best = cuts[np.argmin(chi)]
            left.hist(eig[lam >= best], bins=edges, density=True, histtype="step",
                      color=ACCENT, linewidth=1.8,
                      label=f"SO({2 * n + 1}) excised, best cut drops "
                            f"{drops[np.argmin(chi)]:.0%}")
    right.axvline(arithmetic, color=SERIES, linewidth=1.4,
                  label=f"Gross-Zagier floor: drops {arithmetic:.2%}")
    left.set_xlim(0, s_max)
    left.set_ylim(0, left.get_ylim()[1] * 1.35)
    left.set_xlabel("x of the lowest zero above the central one", color=INK)
    left.set_ylabel("density", color=INK)
    left.set_title("Odd twists of 11a3 at the effective matrix size", color=INK, fontsize=11,
                   loc="left")
    right.set_xlim(0, drops[-1])
    right.set_ylim(0, 250)
    right.set_xlabel("fraction of matrices excised (smallest |Λ'(1)| first)", color=INK)
    right.set_ylabel(f"χ² on {bins} bins of (0, {s_max})", color=INK)
    right.set_title("Excision at N near N_eff", color=INK, fontsize=11, loc="left")
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        ax.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
        _style(ax)
    fig.tight_layout()
    return fig


def _pair_gap_squared(z: np.ndarray, xc: float) -> float:
    upper = z[z.imag > 0]
    if len(upper):
        return float(-4 * upper[np.argmin(np.abs(upper.real - xc))].imag ** 2)
    x = np.sort(z.real)
    i = np.searchsorted(x, xc)
    return float((x[i] - x[i - 1]) ** 2)


def plot_heat_flow(ts: ArrayLike, paths: list[np.ndarray], tc: float, xc: float) -> Figure:
    import matplotlib.pyplot as plt

    ts = np.asarray(ts, dtype=float)
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(12, 4.8), facecolor=SURFACE, gridspec_kw={"width_ratios": [1.35, 1]}
    )
    for t, z in zip(ts, paths, strict=True):
        real = z[z.imag == 0].real
        left.plot(real, np.full(len(real), t), "o", color=SERIES, markersize=2.6)
        up = z[z.imag > 0]
        left.plot(up.real, np.full(len(up), t), "o", color=ACCENT, markersize=2.6)
    left.plot([xc], [tc], "x", color=INK, markersize=8, markeredgewidth=1.6,
              label=f"double zero at t = {tc:.5f}, x = {xc:.4f}")
    left.plot([], [], "o", color=SERIES, markersize=4, label="real zeros")
    left.plot([], [], "o", color=ACCENT, markersize=4, label="Re z of the complex pair")
    left.axhline(0, color=INK_SECONDARY, linewidth=0.8, linestyle=":")
    left.set_ylim(ts.min(), ts.max() + 0.12 * (ts.max() - ts.min()))
    left.set_xlabel("x (zeros of H₀ are x = 2γ)", color=INK)
    left.set_ylabel("t", color=INK)
    left.set_title("Zeros of H_t near the pair γ = 564.161, 564.506", color=INK, fontsize=11,
                   loc="left")
    left.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right",
                ncol=3)
    gap2 = np.array([_pair_gap_squared(z, xc) for z in paths])
    right.plot(ts, gap2, "o", color=SERIES, markersize=3.2, label="(z₊ − z₋)², measured")
    right.plot(ts, 8 * (ts - tc), color=ACCENT, linewidth=1.4, label="8 (t − t_c)")
    right.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    right.axvline(tc, color=INK_SECONDARY, linewidth=0.8, linestyle=":")
    right.set_xlabel("t", color=INK)
    right.set_ylabel("squared gap of the pair", color=INK)
    right.set_title("Real pair, collision, complex pair", color=INK, fontsize=11, loc="left")
    right.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper left")
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def heat_flow_animation(
    ts: ArrayLike, paths: list[np.ndarray], lo: float, hi: float, samples: int = 800,
):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation

    ts = np.asarray(ts, dtype=float)
    x = np.linspace(lo, hi, samples)
    norm = np.exp(np.pi * x / 8) * x ** -1.75
    curves = np.array([np.asarray(ht(x, t)).real * norm for t in ts])
    reach = np.abs(curves).max() * 1.08
    y_max = max(0.1, max(np.abs(z.imag).max() for z in paths)) * 1.3
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(10, 6.2), facecolor=SURFACE, sharex=True,
                                      gridspec_kw={"height_ratios": [1.3, 1]})
    for ax in (top, bottom):
        ax.set_facecolor(SURFACE)
        _style(ax)
    curve, = top.plot([], [], color=SERIES, linewidth=1.4)
    marks, = top.plot([], [], "o", color=ACCENT, markersize=5, markeredgecolor=SURFACE)
    top.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    top.set_xlim(lo, hi)
    top.set_ylim(-reach, reach)
    top.set_ylabel("H_t(x) · e^(πx/8) x^(−7/4)", color=INK)
    dots, = bottom.plot([], [], "o", color=ACCENT, markersize=6, markeredgecolor=SURFACE)
    bottom.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    bottom.set_ylim(-y_max, y_max)
    bottom.set_xlabel("Re z", color=INK)
    bottom.set_ylabel("Im z", color=INK)
    title = fig.suptitle("", color=INK, fontsize=12, x=0.01, ha="left")
    fig.tight_layout()

    def draw(i):
        z = paths[i]
        curve.set_data(x, curves[i])
        real = z[z.imag == 0].real
        marks.set_data(real, np.zeros(len(real)))
        both = np.concatenate([z, np.conj(z[z.imag > 0])])
        dots.set_data(both.real, both.imag)
        state = "all zeros real" if len(real) == len(z) else "a pair has left the real line"
        title.set_text(f"t = {ts[i]:+.3f}: {state}")
        return curve, marks, dots, title

    return FuncAnimation(fig, draw, frames=len(ts), interval=90, blit=False)


def plot_race(
    q: int = 4, a: int = 3, b: int = 1, T: float = 3000.0, x_max: float = 1e6,
    zoom: tuple[float, float] = (26_500.0, 27_200.0), samples: int = 3000,
) -> Figure:
    import matplotlib.pyplot as plt

    R = race(q, a, b, T)
    p, d = prime_race(q, a, b, int(x_max))
    phi = totient(q)
    fig, (left, right) = plt.subplots(
        1, 2, figsize=(12, 4.6), facecolor=SURFACE, gridspec_kw={"width_ratios": [1.6, 1]}
    )
    x = np.exp(np.linspace(np.log(100), np.log(x_max), samples))
    i = np.searchsorted(p, x, side="right") - 1
    left.plot(x, phi * d[i] * np.log(x) / np.sqrt(x), color="#b9b8b3", linewidth=2.2,
              label="primes")
    left.plot(x, R.explicit(x), color=SERIES, linewidth=0.6,
              label=f"{len(R.gamma):,} zeros of L(s, χ)")
    left.axhline(R.bias, color=ACCENT, linewidth=1.2, linestyle="--", label=f"mean {R.bias:g}")
    left.axhline(0, color=INK, linewidth=0.8)
    left.set_xscale("log")
    left.set_xlim(100, x_max)
    left.set_xlabel("x", color=INK)
    left.set_ylabel(f"{phi}·(π(x; {q}, {a}) − π(x; {q}, {b}))·log x / √x", color=INK)
    left.set_title(f"The race between {a} and {b} mod {q}", color=INK, fontsize=11, loc="left")
    xz = np.linspace(*zoom, 1400)
    iz = np.searchsorted(p, xz, side="right") - 1
    right.step(xz, d[iz], where="post", color=INK_SECONDARY, linewidth=1.0, label="primes")
    right.plot(xz, R.explicit(xz) * np.sqrt(xz) / (phi * np.log(xz)), color=SERIES,
               linewidth=1.4, label="zeros")
    right.axhline(0, color=INK, linewidth=0.8)
    first = p[np.argmax(d < 0)]
    if zoom[0] <= first <= zoom[1]:
        right.annotate(f"{first}", (first, -1), xytext=(first + 60, -1.6), color=ACCENT,
                       fontsize=9, arrowprops={"arrowstyle": "-", "color": ACCENT})
    right.set_xlim(*zoom)
    right.set_ylim(-2.5, max(4.5, d[iz].max() + 1))
    right.set_xlabel("x", color=INK)
    right.set_ylabel(f"π(x; {q}, {a}) − π(x; {q}, {b})", color=INK)
    right.set_title(f"First lead for {b} mod {q}", color=INK, fontsize=11, loc="left")
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        ax.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
        _style(ax)
    fig.tight_layout()
    return fig


def plot_bias_densities(races: list[Race], u: ArrayLike | None = None) -> Figure:
    import matplotlib.pyplot as plt

    u = np.linspace(-4, 16, 2001) if u is None else np.asarray(u, dtype=float)
    fig, ax = plt.subplots(figsize=(9, 4.6), facecolor=SURFACE)
    colors = [ORDINAL[1], ORDINAL[3], ACCENT, INK_SECONDARY, SERIES]
    for R, c in zip(races, colors, strict=False):
        f = R.pdf(u)
        ax.plot(u, f, color=c, linewidth=1.6,
                label=f"q = {R.q}, {R.a} vs {R.b}: δ = {R.density():.6f}")
        ax.fill_between(u[u <= 0], f[u <= 0], color=c, alpha=0.25, linewidth=0)
    ax.axvline(0, color=INK, linewidth=0.8)
    ax.set_xlim(u[0], u[-1])
    ax.set_ylim(0, None)
    ax.set_xlabel("limiting value of φ(q)·(π(x; q, a) − π(x; q, b))·log x / √x", color=INK)
    ax.set_ylabel("density", color=INK)
    ax.set_title("Limiting distributions under GRH and LI; shaded where b leads",
                 color=INK, fontsize=11, loc="left")
    ax.set_facecolor(SURFACE)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
    _style(ax)
    fig.tight_layout()
    return fig


def plot_moments(T: ArrayLike, H: float, integrals: ArrayLike) -> Figure:
    import matplotlib.pyplot as plt

    T = np.asarray(T, dtype=float)
    integrals = np.asarray(integrals, dtype=float)
    fig, (full, lead) = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURFACE)
    for k, color in zip((1, 2, 3), (ORDINAL[1], ORDINAL[2], ORDINAL[4]), strict=True):
        data = integrals[:, k - 1]
        cfkrs = np.array([polynomial_integral(k, t, t + H) for t in T])
        mid = np.log((T + H / 2) / (2 * np.pi))
        leading = H * keating_snaith(k, mid)
        for ax, ref in ((full, cfkrs), (lead, leading)):
            ax.plot(T, data / ref, "o-", color=color, markersize=4, linewidth=1.4, label=f"k = {k}")
    full.axhline(1, color=INK_SECONDARY, linewidth=0.8)
    lead.axhline(1, color=INK_SECONDARY, linewidth=0.8)
    lead.set_yscale("log")
    for ax, title in ((full, "against the full CFKRS polynomial"),
                      (lead, "against the leading term only")):
        ax.set_facecolor(SURFACE)
        _style(ax)
        ax.set_xscale("log")
        ax.set_xlabel("T", color=INK)
        ax.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY)
        ax.set_title(title, color=INK, fontsize=12, loc="left")
    full.set_ylabel(f"∫ |ζ|^2k over [T, T + {H:,.0f}] / prediction", color=INK)
    fig.tight_layout()
    return fig


def plot_li(n: ArrayLike, lam: ArrayLike) -> Figure:
    import matplotlib.pyplot as plt

    n = np.asarray(n, dtype=float)
    lam = np.asarray(lam, dtype=float)
    fig, (left, res) = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURFACE)
    left.plot(n, lam, color=SERIES, linewidth=1.6, label="λ_n from the zeros")
    big = n >= 20
    left.plot(n[big], li_trend(n[big]), color=ACCENT, linewidth=1.2, linestyle="--",
             label="(n/2)(log n − 1 + γ − log 2π)")
    left.set_xscale("log")
    left.set_yscale("log")
    left.set_ylabel("λ_n", color=INK)
    left.set_title("Li's coefficients are positive", color=INK, fontsize=12, loc="left")
    res.plot(n, lam - li_trend(n), color=SERIES, linewidth=0.5)
    res.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    res.set_xlim(0, n[-1])
    res.set_ylabel("λ_n − trend", color=INK)
    res.set_title("What is left after the trend", color=INK, fontsize=12, loc="left")
    for ax in (left, res):
        ax.set_facecolor(SURFACE)
        _style(ax)
        ax.set_xlabel("n", color=INK)
    left.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY)
    fig.tight_layout()
    return fig


def _style(ax) -> None:
    ax.grid(color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)


def plot_jensen(
    ns: ArrayLike, roots: ArrayLike, errors: dict[int, ArrayLike], delta: ArrayLike, d: int
) -> Figure:
    import matplotlib.pyplot as plt

    from zeta_explorer.jensen import hermite, real_roots

    ns = np.asarray(ns, dtype=float)
    roots = np.asarray(roots, dtype=float)
    fig, (left, right) = plt.subplots(1, 2, figsize=(12, 4.6), facecolor=SURFACE)
    for x in real_roots(hermite(d)):
        left.axhline(x, color=ACCENT, linewidth=1, linestyle="--")
    for j in range(roots.shape[1]):
        left.plot(ns, roots[:, j], color=SERIES, linewidth=1.6)
    left.set_xscale("log")
    left.set_ylabel("zeros of the renormalized J^(d,n)", color=INK)
    left.set_title(f"d = {d}: zeros move to those of H_{d} (dashed)", color=INK, fontsize=12,
                   loc="left")
    delta = np.asarray(delta, dtype=float)
    colors = (ORDINAL[1], ORDINAL[2], ORDINAL[4])
    for (k, err), color in zip(sorted(errors.items()), colors, strict=False):
        right.plot(ns, np.asarray(err) / delta, "o-", color=color, markersize=3, linewidth=1.4,
                   label=f"d = {k}")
    right.set_xscale("log")
    right.set_ylim(bottom=0)
    right.set_ylabel("largest distance to a zero of H_d, over δ(n)", color=INK)
    right.set_title("The error is a constant times δ(n)", color=INK, fontsize=12, loc="left")
    right.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY)
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
        ax.set_xlabel("shift n", color=INK)
    fig.tight_layout()
    return fig


def plot_maxima(
    zeta: dict[float, ArrayLike], cue: dict[int, ArrayLike], shift: ArrayLike, bins: int = 80
) -> Figure:
    import matplotlib.pyplot as plt

    from zeta_explorer.maxima import fhk_density
    from zeta_explorer.moments import EULER_GAMMA

    fig, (hist, mean, sd) = plt.subplots(1, 3, figsize=(15, 4.6), facecolor=SURFACE)
    shift = np.asarray(shift, dtype=float)
    N_top = max(zeta)
    y = np.asarray(zeta[N_top], dtype=float)
    c = np.asarray(cue[round(N_top)], dtype=float)
    hybrid = c + shift[: len(c)]
    edges = np.linspace(np.quantile(y, 1e-4) - 0.3, np.quantile(y, 1 - 1e-4) + 0.3, bins + 1)
    hist.hist(y, edges, density=True, color=ORDINAL[0], label=f"ζ, N = {N_top:.1f}")
    hist.hist(c, edges, density=True, histtype="step", color=INK_SECONDARY, linewidth=1.2,
              label=f"CUE, N = {round(N_top)}")
    hist.hist(hybrid, edges, density=True, histtype="step", color=ORDINAL[4], linewidth=1.6,
              label="CUE + small primes")
    grid = np.linspace(edges[0], edges[-1], 400)
    hist.plot(grid, fhk_density(grid), color=ACCENT, linewidth=1.6, linestyle="--",
              label="N → ∞ (two Gumbels)")
    hist.set_xlabel("max log|ζ| − log N + ¾ log log N", color=INK)
    hist.set_ylabel("density", color=INK)
    hist.set_title("Maximum over a window of length 2π", color=INK, fontsize=12, loc="left")
    zn, cn = sorted(zeta), sorted(cue)
    for ax, stat, limit, label in ((mean, np.mean, EULER_GAMMA, "mean"),
                                   (sd, np.std, np.pi / np.sqrt(12), "standard deviation")):
        ax.plot(zn, [stat(zeta[n]) for n in zn], "o-", color=SERIES, markersize=4,
                linewidth=1.4, label="ζ")
        ax.plot(cn, [stat(cue[n]) for n in cn], "s-", color=INK_SECONDARY, markersize=4,
                linewidth=1.2, label="CUE")
        ax.axhline(limit, color=ACCENT, linewidth=1.2, linestyle="--", label="N → ∞")
        ax.set_xscale("log")
        ax.set_xlabel("N = log(T/2π), or the matrix size", color=INK)
        ax.set_title(label, color=INK, fontsize=12, loc="left")
    sd.plot(cn, [np.sqrt(np.var(cue[n]) + shift.var()) for n in cn], "s-", color=ORDINAL[4],
            markersize=4, linewidth=1.2, label="CUE + small primes")
    for ax in (hist, mean, sd):
        ax.set_facecolor(SURFACE)
        _style(ax)
    hist.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY, loc="upper right")
    mean.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY)
    sd.legend(frameon=False, fontsize=9, labelcolor=INK_SECONDARY)
    fig.tight_layout()
    return fig


def plot_davenport(off: ArrayLike, line: ArrayLike, sigma0: float, T: float) -> Figure:
    import matplotlib.pyplot as plt

    from zeta_explorer.davenport import dh_z

    off = np.asarray(off, dtype=complex)
    line = np.sort(np.asarray(line, dtype=float))
    fig, (zp, zmap, share) = plt.subplots(1, 3, figsize=(15, 4.8), facecolor=SURFACE)

    t = np.linspace(80, 92, 2401)
    zp.plot(t, dh_z(t), color=SERIES, linewidth=1.5)
    zp.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    near = line[(line > t[0]) & (line < t[-1])]
    zp.plot(near, np.zeros_like(near), "o", color=SERIES, markersize=4)
    pair = off[(off.imag > t[0]) & (off.imag < t[-1]) & (off.real > 0.5)]
    for rho in pair:
        zp.axvline(rho.imag, color=ACCENT, linewidth=1.2, linestyle="--")
        zp.text(0.98, 0.04, f"zeros at {rho.real:.4f} and {1 - rho.real:.4f} + {rho.imag:.3f}i",
                transform=zp.transAxes, ha="right", color=ACCENT, fontsize=9)
    zp.set_xlabel("t", color=INK)
    zp.set_ylabel("Z_f(t)", color=INK)
    zp.set_title("Z_f misses a sign change at the first pair", color=INK, fontsize=12,
                 loc="left")

    for edge in (sigma0, 1 - sigma0):
        zmap.axvline(edge, color=INK_SECONDARY, linewidth=1, linestyle=":")
    zmap.axvspan(sigma0, sigma0 + 0.1, color=GRID)
    zmap.axvspan(1 - sigma0 - 0.1, 1 - sigma0, color=GRID)
    zmap.axvline(0.5, color=SERIES, linewidth=2, alpha=0.35)
    zmap.axvline(1, color=INK_SECONDARY, linewidth=0.8)
    zmap.axvline(0, color=INK_SECONDARY, linewidth=0.8)
    zmap.plot(off.real, off.imag, "o", color=ACCENT, markersize=2.5)
    zmap.set_xlim(1 - sigma0 - 0.1, sigma0 + 0.1)
    zmap.set_ylim(0, T)
    zmap.set_xlabel("σ", color=INK)
    zmap.set_ylabel("t", color=INK)
    zmap.set_title(f"{len(off)} zeros off the line, {len(line)} on it", color=INK, fontsize=12,
                   loc="left")

    grid = np.linspace(200, T, 400)
    n_off = np.searchsorted(np.sort(off.imag), grid)
    n_line = np.searchsorted(line, grid)
    share.plot(grid, n_off / (n_off + n_line), color=ACCENT, linewidth=1.6)
    share.set_ylim(bottom=0)
    share.set_xlim(0, T)
    share.set_xlabel("T", color=INK)
    share.set_ylabel("share of zeros with 0 < t < T off the line", color=INK)
    share.set_title("Share of zeros off the line", color=INK, fontsize=12, loc="left")
    for ax in (zp, zmap, share):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_delta(st: SatoTate, zeros: ArrayLike, x: ArrayLike, ms: tuple[int, ...] = (1, 2, 3),
               bins: int = 60) -> Figure:
    import matplotlib.pyplot as plt

    from zeta_explorer import delta

    zeros = np.asarray(zeros, dtype=float)
    x = np.asarray(x, dtype=float)
    fig, (zp, hist, bias) = plt.subplots(1, 3, figsize=(15, 4.8), facecolor=SURFACE)

    t = np.linspace(0, 40, 1601)
    zp.plot(t, delta.hardy_z(t), color=SERIES, linewidth=1.5)
    zp.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    near = zeros[zeros <= t[-1]]
    zp.plot(near, np.zeros_like(near), "o", color=ACCENT, markersize=4)
    zp.set_xlim(t[0], t[-1])
    zp.set_xlabel("t", color=INK)
    zp.set_ylabel("Z_Δ(t)", color=INK)
    zp.set_title(f"Z for L(s, Δ), first zero {zeros[0]:.6f}", color=INK, fontsize=12, loc="left")

    hist.hist(st.angle, bins=bins, range=(0, np.pi), density=True, color=ORDINAL[0],
              edgecolor=SURFACE, linewidth=0.4)
    th = np.linspace(0, np.pi, 400)
    hist.plot(th, delta.st_density(th), color=ACCENT, linewidth=1.8)
    hist.set_xlim(0, np.pi)
    hist.set_xlabel("θ_p, with τ(p) = 2 p^(11/2) cos θ_p", color=INK)
    hist.set_ylabel("density", color=INK)
    hist.set_title(f"{len(st.p):,} primes against (2/π) sin²θ", color=INK, fontsize=12,
                   loc="left")

    shades = iter(ORDINAL[1::2])
    for m in ms:
        odd = delta.root_number(m) == -1
        bias.plot(x, st.psi(m, x) / np.sqrt(x), linewidth=1.1,
                  color=ACCENT if odd else next(shades),
                  label=f"Sym^{m}, ε = {delta.root_number(m):+d}")
    for level in (0, -2):
        bias.axhline(level, color=INK_SECONDARY, linewidth=1, linestyle="--")
    bias.set_xscale("log")
    bias.set_xlim(x[0], x[-1])
    bias.set_xlabel("x", color=INK)
    bias.set_ylabel("Ψ_m(x) / √x", color=INK)
    bias.set_title("Central zeros of Sym^m L show up in the primes", color=INK, fontsize=12,
                   loc="left")
    bias.set_ylim(-4.5, 2.5)
    bias.legend(frameon=False, fontsize=8, ncol=3, loc="upper center")
    for ax in (zp, hist, bias):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_dedekind(D: int, gz: ArrayLike, gl: ArrayLike, t0: float = 100.0, width: float = 12.0,
                  bins: int = 60) -> Figure:
    import matplotlib.pyplot as plt

    from zeta_explorer import dedekind

    gz, gl = np.asarray(gz, dtype=float), np.asarray(gl, dtype=float)
    chi = quadratic(D)
    fig, (zp, sp) = plt.subplots(1, 2, figsize=(13, 4.8), facecolor=SURFACE)

    t = np.linspace(t0, t0 + width, 1201)
    zp.plot(t, z_rs(t), color=SERIES, linewidth=1.5, label="Z(t) for ζ")
    zp.plot(t, l_z(t, chi), color=ACCENT, linewidth=1.5, label=f"Z(t) for L(s, χ_{D})")
    zp.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    for g, color in ((gz, SERIES), (gl, ACCENT)):
        near = g[(g > t[0]) & (g < t[-1])]
        zp.plot(near, np.zeros_like(near), "o", color=color, markersize=4)
    zp.set_xlim(t[0], t[-1])
    zp.set_xlabel("t", color=INK)
    zp.set_title(f"ζ_K = ζ · L(s, χ_{D}): two independent sets of zeros", color=INK, fontsize=12,
                 loc="left")
    lo, hi = zp.get_ylim()
    zp.set_ylim(lo, hi + 0.35 * (hi - lo))
    zp.legend(frameon=False, fontsize=9, ncol=2, loc="upper left")

    u = np.sort(np.concatenate([dedekind.unfold(D, gz), dedekind.unfold(D, gl)]))
    s = np.diff(u)
    sp.hist(s, bins=bins, range=(0, 3), density=True, color=ORDINAL[0], edgecolor=SURFACE,
            linewidth=0.4)
    x = np.linspace(0, 3, 400)
    sp.plot(x, dedekind.superposition_spacing(x), color=ACCENT, linewidth=1.8,
            label="two independent GUE halves")
    sp.plot(x, gue_spacing(x), color=INK_SECONDARY, linewidth=1.4, linestyle="--",
            label="one GUE spectrum")
    sp.set_xlim(0, 3)
    sp.set_xlabel("normalized spacing", color=INK)
    sp.set_ylabel("density", color=INK)
    sp.set_title(f"{len(s):,} spacings of ζ_K zeros up to t = {max(gz[-1], gl[-1]):,.0f}",
                 color=INK, fontsize=12, loc="left")
    sp.legend(frameon=False, fontsize=9)
    for ax in (zp, sp):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_mertens(x: ArrayLike, q: ArrayLike, qz: ArrayLike, t: ArrayLike,
                 profiles: list[tuple[str, np.ndarray]], ns: ArrayLike, h_pos: ArrayLike,
                 h_neg: ArrayLike, zeros: int) -> Figure:
    import matplotlib.pyplot as plt

    x, q, qz, t = (np.asarray(v, dtype=float) for v in (x, q, qz, t))
    fig = plt.figure(figsize=(13, 8.4), facecolor=SURFACE)
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.05])
    top = fig.add_subplot(grid[0, :])
    left, right = fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])

    top.plot(x, q, color=SERIES, linewidth=2.2, label="M(x)/√x from the sieve")
    top.plot(x, qz, color=ACCENT, linewidth=0.6, label=f"explicit formula, {zeros:,} zeros")
    for v in (1, -1):
        top.axhline(v, color=INK_SECONDARY, linewidth=0.8, linestyle="--")
    top.set_xscale("log")
    top.set_xlim(x[0], x[-1])
    top.set_ylim(-1.15, 1.45)
    top.set_xlabel("x", color=INK)
    top.set_title("Mertens: |M(x)| < √x as far as anyone can sieve", color=INK, fontsize=12,
                  loc="left")
    top.legend(frameon=False, fontsize=9, ncol=2, loc="upper left")

    for (label, p), color in zip(profiles, (ACCENT, SERIES), strict=False):
        left.plot(t, p, color=color, linewidth=1.1, label=label)
    for v in (1, -1):
        left.axhline(v, color=INK_SECONDARY, linewidth=0.8, linestyle="--")
    left.set_xlim(t[0], t[-1])
    left.set_xlabel("log x − y±", color=INK)
    left.set_ylabel("h(log x)", color=INK)
    left.set_title("Near the y found by LLL the smoothed sum passes ±1", color=INK, fontsize=12,
                   loc="left")
    lo, hi = left.get_ylim()
    left.set_ylim(lo, hi + 0.3 * (hi - lo))
    left.legend(frameon=False, fontsize=9, loc="upper left")

    ns = np.asarray(ns, dtype=float)
    right.plot(ns, h_pos, "o-", color=ACCENT, linewidth=1.4, markersize=4,
               label="lim sup M(x)/√x ≥ h(y₊)")
    right.plot(ns, -np.asarray(h_neg, dtype=float), "s-", color=SERIES, linewidth=1.4,
               markersize=4, label="lim inf M(x)/√x ≤ h(y₋), sign flipped")
    right.axhline(1, color=INK_SECONDARY, linewidth=0.8, linestyle="--")
    right.set_xlabel("zeros in the lattice", color=INK)
    right.set_ylabel("|h(y)|", color=INK)
    right.set_title("Bigger lattices, bigger values", color=INK, fontsize=12, loc="left")
    right.legend(frameon=False, fontsize=9, loc="lower right")
    for ax in (top, left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_skewes(omega: ArrayLike, s: ArrayLike, peaks: ArrayLike, w: ArrayLike, H: ArrayLike,
                bar: float, best: float, zeros: int) -> Figure:
    import matplotlib.pyplot as plt

    omega, s, w, H = (np.asarray(v, dtype=float) for v in (omega, s, w, H))
    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 4.8), facecolor=SURFACE)
    left.plot(omega, s, color=SERIES, linewidth=0.6)
    left.axhline(1, color=INK_SECONDARY, linewidth=0.8, linestyle="--")
    for k, p in enumerate(peaks):
        i = int(np.argmin(np.abs(omega - p)))
        left.plot(omega[i], s[i], "o", color=ACCENT, markersize=5)
        left.annotate(f"{p:.3f}", (omega[i], s[i]), textcoords="offset points",
                      xytext=(-6 if k % 2 == 0 else 6, 14), ha="right" if k % 2 == 0 else "left",
                      fontsize=9, color=ACCENT)
    left.set_xlim(omega[0], omega[-1])
    left.set_ylim(top=1.16)
    left.set_xlabel("ω = log x", color=INK)
    left.set_ylabel("S(ω)", color=INK)
    left.set_title("Lehman's sum over the first 1,000 zeros, max per unit of ω", color=INK,
                   fontsize=12, loc="left")

    base = round(best, 2)
    right.plot(w - base, H, color=SERIES, linewidth=1.2, label=f"H(ω) with {zeros:,} zeros")
    right.axhline(bar, color=ACCENT, linewidth=1.2, linestyle="--", label="1 + |R| + |H − H*|")
    right.fill_between(w - base, bar, H, where=H > bar, color=ACCENT, alpha=0.25, linewidth=0)
    right.axvline(best - base, color=INK_SECONDARY, linewidth=0.8)
    right.set_xlim(w[0] - base, w[-1] - base)
    right.set_xlabel(f"ω − {base:g}", color=INK)
    right.set_title("Above the dashed line, π(x) > li(x) somewhere near $e^ω$", color=INK,
                    fontsize=12, loc="left")
    right.legend(frameon=False, fontsize=9, loc="upper left")
    for ax in (left, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_effective(x: ArrayLike, err: dict[str, ArrayLike], bound: dict[str, ArrayLike],
                   t: float, y: float, rms: dict[float, tuple[ArrayLike, ArrayLike]],
                   model: dict[float, tuple[ArrayLike, ArrayLike]],
                   hist: dict[float, tuple[ArrayLike, ArrayLike]], X: float,
                   bins: int = 50) -> Figure:
    import matplotlib.pyplot as plt

    fig, (left, mid, right) = plt.subplots(1, 3, figsize=(15, 4.6), facecolor=SURFACE)
    for (name, e), color in zip(err.items(), (ORDINAL[0], SERIES), strict=True):
        left.plot(x, e, "o", color=color, markersize=3.5, label=f"|H/B − ({name})|")
        left.plot(x, bound[name], color=color, linewidth=1.4, label=f"bound for {name}")
    left.set_xscale("log")
    left.set_yscale("log")
    left.set_xlabel("x", color=INK)
    left.set_ylabel("error relative to B_t", color=INK)
    left.set_title(f"Effective approximation, t = {t:g}, y = {y:g}", color=INK, fontsize=12,
                   loc="left")
    left.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="lower left")
    colors = dict(zip(sorted(rms), ["#9a9893", *ORDINAL[1:]], strict=False))
    for ts in sorted(rms):
        h, r = rms[ts]
        mid.plot(h, r, "o", color=colors[ts], markersize=4.5, label=f"t = {ts:g}")
        if ts in model:
            g, m = model[ts]
            mid.plot(g, m, color=colors[ts], linewidth=1.3)
    mid.plot([], [], color=INK_SECONDARY, linewidth=1.3, label="random phase model")
    mid.set_xscale("log")
    mid.set_ylim(0, None)
    mid.set_xlabel("x", color=INK)
    mid.set_ylabel("rms displacement from the ideal point", color=INK)
    mid.set_title("Zeros of H_t crystallize", color=INK, fontsize=12, loc="left")
    mid.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="lower left", ncol=2)
    edges = np.linspace(-0.5, 0.5, bins + 1)
    for ts in sorted(hist):
        d, m = hist[ts]
        right.hist(d, edges, density=True, histtype="step", color=colors[ts], linewidth=1.6,
                   label=f"t = {ts:g}")
        if m is not None:
            dens, _ = np.histogram(m, edges, density=True)
            right.plot((edges[1:] + edges[:-1]) / 2, dens, "--", color=colors[ts],
                       linewidth=1.0)
    right.plot([], [], "--", color=INK_SECONDARY, linewidth=1.0, label="model")
    right.set_xlabel("displacement (in units of the mean gap)", color=INK)
    right.set_ylabel("density", color=INK)
    right.set_title(f"Displacements near x = {_fmt_power(X)}", color=INK, fontsize=12,
                    loc="left")
    right.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
    for ax in (left, mid, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_sympower(
    curves: dict[str, tuple[ArrayLike, ArrayLike, ArrayLike]], t: ArrayLike, phi: ArrayLike,
    zeros: ArrayLike, window: tuple[float, float], name: str,
) -> Figure:
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(12, 6.4), facecolor=SURFACE)
    grid = fig.add_gridspec(len(curves), 2, width_ratios=[1.25, 1])
    left = [fig.add_subplot(grid[i, 0]) for i in range(len(curves))]
    for ax, (label, (ts, z, g)), color in zip(left, curves.items(), ORDINAL, strict=False):
        ax.axhline(0, color=INK_SECONDARY, linewidth=0.6)
        ax.plot(ts, z, color=color, linewidth=1.3)
        ax.plot(g, np.zeros(len(g)), "o", color=INK, markersize=2.8)
        ax.set_xlim(ts[0], ts[-1])
        ax.set_ylabel(label, color=INK)
        ax.set_facecolor(SURFACE)
        _style(ax)
    for ax in left[:-1]:
        ax.tick_params(labelbottom=False)
    left[0].set_title("Z(t) for symmetric powers of Δ, zeros marked", color=INK, fontsize=11,
                      loc="left")
    left[-1].set_xlabel("t", color=INK)
    right = fig.add_subplot(grid[:, 1])
    t = np.asarray(t, dtype=float)
    g = np.sort(np.asarray(zeros, dtype=float))
    s = np.searchsorted(g, t, side="right") - np.asarray(phi, dtype=float)
    right.axvspan(*window, color=ACCENT, alpha=0.12, label="Turing window")
    right.axhline(0, color=INK_SECONDARY, linewidth=0.6)
    right.plot(t, s, color=SERIES, linewidth=1.1, label="zeros found − Φ(t)")
    right.plot(t, s - (t > g[len(g) // 3]), color=INK_SECONDARY, linewidth=0.9, linestyle=":",
               label="the same with one zero left out")
    right.set_xlim(t[0], t[-1])
    right.set_ylim(-2.4, 2.2)
    right.set_xlabel("t", color=INK)
    right.set_title(f"S(t) for {name}", color=INK, fontsize=11, loc="left")
    right.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper left")
    right.set_facecolor(SURFACE)
    _style(right)
    fig.tight_layout()
    return fig


def plot_barrier(images: dict[float, ArrayLike], ts: ArrayLike, dts: ArrayLike, points: ArrayLike,
                 edges: list[tuple[float, float, dict[str, float]]], threshold: float,
                 X: float) -> Figure:
    import matplotlib.pyplot as plt

    fig, (left, mid, right) = plt.subplots(1, 3, figsize=(15, 4.8), facecolor=SURFACE)
    for (t, f), color in zip(sorted(images.items()), ORDINAL[::2], strict=False):
        f = np.asarray(f)
        left.plot(np.append(f.real, f.real[0]), np.append(f.imag, f.imag[0]), color=color,
                  linewidth=1.3, label=f"t = {t:g}")
    left.plot([0], [0], "+", color=ACCENT, markersize=11, markeredgewidth=1.8, label="origin")
    left.set_aspect("equal", adjustable="datalim")
    left.set_xlabel("Re f_t", color=INK)
    left.set_ylabel("Im f_t", color=INK)
    left.set_title(f"f_t on the barrier rectangle at x = {_fmt_power(X)}", color=INK,
                   fontsize=12, loc="left")
    left.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper left")
    ts = np.asarray(ts, dtype=float)
    mid.stairs(points, np.append(ts, ts[-1] + np.asarray(dts)[-1]), color=SERIES, linewidth=1.4)
    mid.set_yscale("log")
    mid.set_xlim(0, None)
    mid.set_xlabel("t", color=INK)
    mid.set_ylabel("mesh points on the rectangle", color=INK)
    mid.set_title(f"{len(ts)} time steps, winding number 0 at each", color=INK,
                  fontsize=12, loc="left")
    styles = {"Lemma 8.5, N tracked": (SERIES, "-"), "Lemma 8.5, coefficients at N+": (
        ORDINAL[0], "--"), "triangle inequality": (INK_SECONDARY, "-")}
    for name, (color, ls) in styles.items():
        for i, (lo, hi, vals) in enumerate(edges):
            right.plot([lo, hi], [vals[name]] * 2, color=color, linestyle=ls, linewidth=2.0,
                       label=name if i == 0 else None)
    right.axhline(threshold, color=ACCENT, linewidth=1.0, linestyle=":",
                  label=f"needed: {threshold:g}")
    right.axhline(0, color=INK_SECONDARY, linewidth=0.6)
    right.set_xscale("log")
    right.set_xlabel("N on the lower edge y = 0.2", color=INK)
    right.set_ylabel("lower bound on dist(E f, (−∞, 0])", color=INK)
    right.set_title("Euler-mollified lower bound", color=INK, fontsize=12, loc="left")
    right.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="center right")
    for ax in (left, mid, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig


def plot_funcfield(lowest: ArrayLike, so: ArrayLike, excised: ArrayLike, k: ArrayLike,
                   scan: list, q: int, N: int) -> Figure:
    import matplotlib.pyplot as plt

    fig, (left, mid, right) = plt.subplots(1, 3, figsize=(15, 4.8), facecolor=SURFACE)
    for x, color, name in [(so, INK_SECONDARY, f"SO({N}), Haar"),
                           (excised, ORDINAL[1], f"SO({N}), det(I − A) ≥ 4/{q}^{N // 4}"),
                           (lowest, ACCENT, f"{len(lowest)} rank 0 twists")]:
        x = np.sort(np.asarray(x))
        left.plot(x, np.arange(1, len(x) + 1) / len(x), color=color, linewidth=1.6, label=name,
                  drawstyle="steps-post" if color == ACCENT else "default")
    left.set_xlim(0, 1.2)
    left.set_ylim(0, 1)
    left.set_xlabel("lowest eigenangle θ₁ (radians, no unfolding)", color=INK)
    left.set_ylabel("cumulative fraction", color=INK)
    left.set_title(f"Even twists of the Legendre curve over F_{q}(t)", color=INK, fontsize=12,
                   loc="left")
    left.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="lower right")
    ks, counts = np.unique(np.asarray(k), return_counts=True)
    mid.bar(ks, counts, width=0.7, color=[ACCENT if v == 0 else SERIES for v in ks])
    mid.set_xlabel(f"k, where {q}^{N // 2} · L(E_D, 1/{q}) = 4 · {q}^{N // 4} · k²", color=INK)
    mid.set_ylabel("even twists", color=INK)
    mid.set_title("Central values are 4 q^(N/4) times a square", color=INK, fontsize=12,
                  loc="left")
    mid.set_xticks(ks)
    qs = [f.q for f in scan]
    right.plot(qs, [f.ks_so for f in scan], "o-", color=INK_SECONDARY, linewidth=1.4,
               label=f"against SO({scan[0].N})")
    right.plot(qs, [f.ks_excised for f in scan], "o-", color=ACCENT, linewidth=1.4,
               label=f"against SO({scan[0].N}) with det(I − A) ≥ 4/q")
    right.plot(qs, [1.36 / np.sqrt(f.even) for f in scan], ":", color=SERIES, linewidth=1.2,
               label="5% critical value")
    right.set_ylim(0, None)
    right.set_xlabel("q", color=INK)
    right.set_ylabel("KS distance, lowest eigenangle", color=INK)
    right.set_title(f"Degree {scan[0].N}: rank 0 even twists as q grows", color=INK,
                    fontsize=12, loc="left")
    right.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right")
    for ax in (left, mid, right):
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig



def plot_symfamily(F: SymFamily, buckets: list[tuple[int, int]],
                   smooth: dict[str, tuple[np.ndarray, ...]], s0: float,
                   sweep: tuple[np.ndarray, ...], x_max: float = 3.0, bins: int = 30) -> Figure:
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(13, 9.4), facecolor=SURFACE)
    (dens, low), (smo, sup) = axes
    x = np.linspace(0, x_max, 400)
    for harm, color, name in ((False, SERIES, "natural average"),
                              (True, ACCENT, "harmonic (Petersson) weights")):
        c, d = F.density(bins, x_max, harmonic=harm)
        dens.step(c, d, where="mid", color=color, linewidth=1.8, label=name)
    dens.plot(x, w_symplectic(x), color=INK, linewidth=1.4, label="USp: 1 − sin 2πx / 2πx")
    dens.plot(x, w_unitary(x), color=INK_SECONDARY, linewidth=1.2, linestyle=":", label="U: 1")
    dens.plot(x, 1 + np.sinc(2 * x), color=INK_SECONDARY, linewidth=1.2, linestyle="--",
              label="SO(even): 1 + sin 2πx / 2πx")
    dens.set_xlim(0, x_max)
    dens.set_ylim(0, 2.2)
    dens.set_xlabel("x (zero height in units of the mean spacing)", color=INK)
    dens.set_ylabel("zeros per L-function per unit x", color=INK)
    dens.set_title(f"L(s, sym² f): {F.size:,} level 1 eigenforms, weight 12 to {F.k.max()}",
                   color=INK, fontsize=11, loc="left")
    dens.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="upper right",
                ncol=2)

    s = np.linspace(0, 1.6, 301)
    for (lo, hi), color in zip(buckets, ORDINAL[::2], strict=False):
        first = np.sort(F.select((F.k >= lo) & (F.k <= hi)).first())
        low.step(first, 1 - np.arange(len(first)) / len(first), where="post", color=color,
                 linewidth=1.8, label=f"k = {lo} to {hi} ({len(first):,} forms)")
    for n, style in ((2, "--"), (6, "-.")):
        low.plot(s, gap(symplectic_n(n), s), color=INK_SECONDARY, linewidth=1.1, linestyle=style,
                 label=f"USp({2 * n})")
    low.plot(s, gap(symplectic, s), color=INK, linewidth=1.4, label="USp, N → ∞")
    low.set_yscale("log")
    low.set_xlim(0, 1.6)
    low.set_ylim(1e-3, 1.2)
    low.set_xlabel("s", color=INK)
    low.set_ylabel("P(lowest zero > s)", color=INK)
    low.set_title("Lowest zero: narrower than USp(2N) for every N", color=INK, fontsize=11,
                  loc="left")
    low.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="lower left")

    for (label, (x0, data, pred, usp)), color in zip(smooth.items(), ORDINAL[::2], strict=False):
        smo.plot(x0, data, "o", color=color, markersize=4, label=f"{label}: zeros")
        smo.plot(x0, pred, color=color, linewidth=1.6, label="explicit formula, E λ(n) = δ(n, 1)")
        smo.plot(x0, usp, color=color, linewidth=1.2, linestyle=":", label="USp limit")
    smo.set_xlim(0, x0.max())
    smo.set_xlabel("x₀", color=INK)
    smo.set_ylabel("Σ φ(x_γ − x₀) + φ(x_γ + x₀), harmonic average", color=INK)
    smo.set_title(f"Smoothed density, Fourier support {s0:g}", color=INK, fontsize=11, loc="left")
    smo.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="lower right")

    supports, harm, nat, lim = sweep
    sup.plot(supports, harm, "o-", color=ACCENT, linewidth=1.6,
             label="harmonic, against E λ(n) = δ(n, 1)")
    sup.plot(supports, nat, "s-", color=SERIES, linewidth=1.6,
             label="natural, against Plancherel E λ(p²ʲ) = p⁻ʲ")
    sup.plot(supports, lim, "^:", color=INK_SECONDARY, linewidth=1.2, label="harmonic, against USp")
    sup.axvline(1.5, color=GRID, linewidth=6, zorder=0)
    sup.set_yscale("log")
    sup.set_xlabel("support of the Fourier transform of φ", color=INK)
    sup.set_ylabel("max |zeros − prediction| over x₀ ∈ [0, 3]", color=INK)
    sup.set_title("Where the diagonal stops (ILS prove 3/2)", color=INK, fontsize=11, loc="left")
    sup.legend(frameon=False, fontsize=8.5, labelcolor=INK_SECONDARY, loc="lower right")
    for ax in axes.ravel():
        ax.set_facecolor(SURFACE)
        _style(ax)
    fig.tight_layout()
    return fig
