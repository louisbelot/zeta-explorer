import math

import matplotlib
import numpy as np
import pytest
from matplotlib.animation import PillowWriter
from PIL import Image

from zeta_explorer.bias import race
from zeta_explorer.dedekind import zeros as field_zeros
from zeta_explorer.delta import sato_tate, zeros
from zeta_explorer.li import li_trend
from zeta_explorer.lowlying import OneLevel, w_symplectic
from zeta_explorer.moments import EULER_GAMMA, polynomial_integral
from zeta_explorer.plot import (
    heat_flow_animation,
    phase_colors,
    plot_afe,
    plot_bias_densities,
    plot_central_values,
    plot_continuation,
    plot_davenport,
    plot_dedekind,
    plot_delta,
    plot_domain,
    plot_ec_z,
    plot_explicit,
    plot_first_zero,
    plot_gram,
    plot_heat_flow,
    plot_jensen,
    plot_l_z,
    plot_l_zeros,
    plot_landscape,
    plot_lehmer,
    plot_li,
    plot_log_z,
    plot_low_lying,
    plot_low_lying_scales,
    plot_maxima,
    plot_moments,
    plot_multieval,
    plot_odd_excision,
    plot_race,
    plot_rs_error,
    plot_spacing,
    plot_spiral,
    plot_twist_zeros,
    plot_z,
    sweep_sigma,
)
from zeta_explorer.rmt import so_angles
from zeta_explorer.spacing import first_zeros, zeros_between
from zeta_explorer.values import log_z_sample

matplotlib.use("Agg")


def test_plot_z_marks_every_zero(tmp_path):
    fig = plot_z(10, 30, samples=200)
    zero_markers = fig.axes[0].lines[-1]
    assert len(zero_markers.get_xdata()) == 3
    fig.savefig(tmp_path / "z.png")
    assert (tmp_path / "z.png").stat().st_size > 0


def test_plot_gram_flags_g126_and_g134():
    fig = plot_gram(262, 300, samples=200)
    (bad,) = [line for line in fig.axes[0].lines if line.get_label() == "bad Gram point"]
    assert len(bad.get_xdata()) == 2


def test_plot_rs_error_draws_one_line_per_order():
    fig = plot_rs_error(200, 400, bins=2)
    assert len(fig.axes[0].lines) == 5


def test_phase_colors_hue_and_bands():
    rgb = phase_colors(np.array([1.0, 1j, -1.0, 2.0, 0.0, np.inf]))
    assert rgb.shape == (6, 3)
    assert rgb[0].argmax() == 0
    assert rgb[2].argmax() != 0
    assert np.allclose(rgb[3], rgb[0])
    assert np.isfinite(rgb).all()


def test_plot_domain_image_size():
    fig = plot_domain(0, 10, -1, 2, width=100)
    assert fig.axes[0].images[0].get_array().shape == (30, 100, 3)


def test_plot_spiral_counts_zeros():
    fig = plot_spiral(t_max=26, samples=300)
    assert any("3 times" in text.get_text() for text in fig.axes[0].texts)


def test_plot_continuation_draws_three_series_twice():
    fig = plot_continuation(N=50)
    left, right = fig.axes
    assert len(right.lines) == 3
    assert len(left.lines) == 4


def test_plot_spacing():
    fig = plot_spacing(zeros_between(1000, 1400), bins=20)
    assert len(fig.axes) == 2
    assert "399 zero spacings" in fig.axes[0].get_legend().get_texts()[0].get_text()


def test_plot_explicit_draws_each_truncation():
    fig = plot_explicit(first_zeros(20), counts=(5, 10, 20), x_max=30, samples=100)
    left, right = fig.axes
    labels = [line.get_label() for line in left.lines]
    assert labels[-3:] == ["5 zeros", "10 zeros", "20 zeros"]
    assert len(right.lines) == 4


def test_plot_lehmer_marks_the_pair():
    fig = plot_lehmer(samples=200)
    zoom = fig.axes[1]
    assert "0.0377 apart" in zoom.get_title(loc="left")
    assert len(zoom.lines[-1].get_xdata()) == 2


def test_plot_log_z():
    T = 2 * math.pi * math.exp(8)
    fig = plot_log_z(log_z_sample(T, 2000, np.random.default_rng(1)), T, 8, bins=30)
    texts = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert texts[0].startswith("log|ζ(1/2 + it)|, 2,000 random t")
    assert "CUE, N = 8" in texts[1]
    assert "arithmetic factor, variance 1.741" in texts[2]
    assert len(fig.axes) == 2


def test_plot_moments_ratios():
    T = np.array([1e4, 1e5])
    rows = np.array([[polynomial_integral(k, t, t + 100) for k in (1, 2, 3)] for t in T])
    fig = plot_moments(T, 100.0, rows)
    for line in fig.axes[0].lines[:3]:
        assert np.allclose(line.get_ydata(), 1)
    assert fig.axes[1].get_yscale() == "log"
    assert fig.axes[1].lines[2].get_ydata()[0] > 10


def test_plot_li_residual():
    n = np.arange(1, 50)
    fig = plot_li(n, li_trend(n) + 1)
    assert np.allclose(fig.axes[1].lines[0].get_ydata(), 1)


def test_plot_landscape_puts_a_pit_at_each_zero():
    fig = plot_landscape(t_max=30, width=120)
    pits = fig.axes[0].collections[-1]
    assert np.allclose(sorted(pits._offsets3d[0]), [14.134725, 21.022040, 25.010858], atol=1e-6)


def test_sweep_sigma_writes_every_frame(tmp_path):
    anim = sweep_sigma(t_max=30, frames=8, samples=300)
    anim.save(tmp_path / "sweep.gif", writer=PillowWriter(fps=6), dpi=40)
    assert Image.open(tmp_path / "sweep.gif").n_frames == 8


def test_plot_multieval_window_beats_direct_sum():
    fig = plot_multieval(heights=(1e4, 1e8), points=2)
    left = fig.axes[0]
    direct, window = (line.get_ydata() for line in left.lines)
    assert np.all(window < 1e-11)
    assert direct[-1] > 1e-9


def test_plot_afe_labels_term_counts():
    fig = plot_afe(conductors=(101, 10007), heights=(1, 100), points=2)
    left, right = fig.axes
    assert [t.get_text() for t in left.texts] == ["100 terms", "986 terms"]
    em, smooth = (line.get_ydata() for line in right.lines)
    assert smooth[-1] < em[-1]


def test_plot_l_z_marks_lmfdb_zeros():
    fig = plot_l_z(labels=((4, 3),), t_min=0, t_max=13, samples=200)
    zeros = fig.axes[0].lines[-1].get_xdata()
    assert np.allclose(zeros, [6.020948904697597, 10.243770304166555, 12.988098012312423])


def test_plot_l_zeros_puts_conjugates_side_by_side():
    fig = plot_l_zeros(q=7, T=10)
    labels = [t.get_text() for t in fig.axes[0].get_xticklabels()]
    assert labels == ["2", "4", "3", "5", "6"]


def test_plot_low_lying(tmp_path):
    a = OneLevel([np.array([0.6, 1.4]), np.array([0.9])])
    b = OneLevel([np.array([0.1, 1.2]), np.array([0.4, 2.2])])
    fig = plot_low_lying(a, b, ("quadratic", "mod p"), bins=6)
    fig.savefig(tmp_path / "low.png")
    assert len(fig.axes) == 2


def test_plot_low_lying_scales_fits_through_origin():
    fams = {1e3: OneLevel([np.array([0.1, 0.5])] * 4), 1e5: OneLevel([np.array([0.2, 0.9])] * 4)}
    fig = plot_low_lying_scales(fams, window=(0.0, 0.3), bins=6)
    excess = 1 / 0.3 - np.mean(w_symplectic(np.linspace(0, 0.3, 2001)))
    u = 1 / np.log([1e3, 1e5])
    slope = excess * u.sum() / np.square(u).sum()
    assert fig.axes[1].get_legend().get_texts()[0].get_text() == f"{slope:.2f} / log |d|"
    assert len(fig.axes[0].lines) == 3


def test_plot_low_lying_scales_draws_ratios():
    fams = {4e3: OneLevel([np.array([0.2, 1.5])] * 4), 1e5: OneLevel([np.array([0.25, 1.6])] * 4)}
    ds = {4e3: np.array([-4003, 4001]), 1e5: np.array([-100003, 100001])}
    fig = plot_low_lying_scales(fams, bins=6, ratios=ds)
    labels = [t.get_text() for t in fig.axes[1].get_legend().get_texts()]
    assert labels[1:] == ["ratios conjecture, d = ±D", "ratios conjecture, same characters"]
    assert len(fig.axes[0].lines) == 6


def test_plot_first_zero(tmp_path):
    fams = {4e3: OneLevel([np.array([0.5]), np.array([0.7, 1.2]), np.array([0.3])])}
    fig = plot_first_zero(fams, n_grid=[1.5, 2.0, 3.0])
    fig.savefig(tmp_path / "first.png")
    assert "N = " in fig.axes[1].get_legend().get_texts()[0].get_text()


def test_plot_ec_z(tmp_path):
    fig = plot_ec_z(t_max=3.0, samples=61)
    fig.savefig(tmp_path / "ec.png")
    titles = [ax.get_title(loc="left") for ax in fig.axes]
    assert [t.split(":")[0] for t in titles] == ["11a3", "37a1", "389a1", "5077a1"]
    assert [t.split("order ")[1][0] for t in titles] == ["0", "1", "2", "3"]


def test_plot_central_values(tmp_path):
    fig = plot_central_values([5, 12, 232], [2.838, 1.832, 0.0], 6.346046521, c_max=3)
    fig.savefig(tmp_path / "central.png")
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert labels == ["L(E_d, 1/2) ≠ 0 (2)", "L(E_d, 1/2) = 0, rank 2 (1)"]


def test_plot_twist_zeros(tmp_path):
    rng = np.random.default_rng(0)
    fig = plot_twist_zeros(rng.uniform(0.1, 1.0, 50), rng.uniform(0.3, 1.4, 50),
                           rng.uniform(0.1, 1.0, 200), 8, 0.0389, n_grid=[1.5, 2.0, 3.0])
    fig.savefig(tmp_path / "twists.png")
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert labels[-1] == "SO(16) with Λ(1) ≥ 0.0389"


def test_plot_odd_excision_finds_the_repulsive_cut(tmp_path):
    theta = so_angles(3, 20_000, seed=1, odd=True)
    lam = np.prod(2 - 2 * np.cos(theta), axis=1) * np.exp(1.5)
    x = 3 * theta[lam >= 20, 0][:800] / np.pi
    fig = plot_odd_excision(x, theta, {"densities": 0.3, "probabilities": 2.0}, c_max=200)
    fig.savefig(tmp_path / "odd.png")
    best = fig.axes[1].get_legend().get_texts()[-1].get_text()
    assert 10 <= float(best.split("= ")[1]) <= 40


def test_plot_heat_flow_squares_the_gap(tmp_path):
    ts = np.array([0.1, 0.0, -0.1])
    paths = [np.array([1.0, 2.0, 4.0 + 0j]), np.array([1.5, 2.5, 4.0 + 0j]),
             np.array([4.0 + 0j, 2.0 + 0.5j])]
    fig = plot_heat_flow(ts, paths, -0.05, 2.0)
    fig.savefig(tmp_path / "heat.png")
    assert list(fig.axes[1].lines[0].get_ydata()) == pytest.approx([1.0, 1.0, -1.0])


def test_heat_flow_animation_writes_every_frame(tmp_path):
    ts = np.array([0.05, 0.0])
    paths = [np.array([28.2 + 0j, 42.1 + 0j]), np.array([28.27 + 0j, 42.05 + 0j])]
    anim = heat_flow_animation(ts, paths, 20.0, 45.0, samples=100)
    anim.save(tmp_path / "heat.gif", writer=PillowWriter(fps=5))
    assert Image.open(tmp_path / "heat.gif").n_frames == 2


def test_plot_race_marks_first_lead(tmp_path):
    fig = plot_race(T=200.0, x_max=30_000, samples=300)
    fig.savefig(tmp_path / "race.png")
    texts = [t.get_text() for t in fig.axes[1].texts]
    assert texts == ["26861"]


def test_plot_bias_densities(tmp_path):
    fig = plot_bias_densities([race(3, 2, 1, 50.0), race(4, 3, 1, 50.0)])
    fig.savefig(tmp_path / "bias.png")
    labels = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert [x.split(":")[0] for x in labels] == ["q = 3, 2 vs 1", "q = 4, 3 vs 1"]


def test_plot_jensen_ratio_panel():
    ns = np.array([1.0, 10.0, 100.0])
    roots = np.tile([-1.0, 1.0], (3, 1))
    delta = np.array([0.4, 0.2, 0.1])
    fig = plot_jensen(ns, roots, {2: 3 * delta}, delta, 2)
    assert np.allclose(fig.axes[1].lines[0].get_ydata(), 3)
    assert len(fig.axes[0].lines) == 4


def test_plot_maxima_limits():
    rng = np.random.default_rng(0)
    zeta = {21.3: rng.gumbel(size=500) / 2 + 0.3}
    cue = {21: rng.gumbel(size=500) / 2, 50: rng.gumbel(size=500) / 2}
    shift = np.full(500, 0.2)
    shift[::2] = -0.2
    fig = plot_maxima(zeta, cue, shift)
    assert fig.axes[1].lines[-1].get_ydata()[0] == EULER_GAMMA
    assert fig.axes[2].lines[2].get_ydata()[0] == pytest.approx(np.pi / np.sqrt(12))
    added = fig.axes[2].lines[3].get_ydata() ** 2 - fig.axes[2].lines[1].get_ydata() ** 2
    assert np.allclose(added, 0.04)


def test_plot_davenport_counts():
    off = np.array([0.8 + 85.7j, 0.2 + 85.7j, 0.6 + 300j, 0.4 + 300j])
    line = np.linspace(5, 400, 300)
    fig = plot_davenport(off, line, 1.4, 400.0)
    assert fig.axes[1].get_title(loc="left") == "4 zeros off the line, 300 on it"
    share = fig.axes[2].lines[0].get_ydata()
    assert np.all((share > 0) & (share < 0.05))
    assert "0.8000" in fig.axes[0].texts[0].get_text()


def test_plot_delta_panels():
    st = sato_tate(1 << 12)
    fig = plot_delta(st, zeros(20), np.geomspace(100, 4000, 50))
    assert fig.axes[0].get_title(loc="left").endswith("9.222379")
    assert fig.axes[1].get_title(loc="left").startswith(f"{len(st.p):,} primes")
    labels = [t.get_text() for t in fig.axes[2].get_legend().get_texts()]
    assert labels == ["Sym^1, ε = +1", "Sym^2, ε = +1", "Sym^3, ε = -1"]


def test_plot_dedekind_counts():
    gz, gl = field_zeros(-4, 120.0, per_gap=16)
    fig = plot_dedekind(-4, gz, gl)
    n = len(gz) + len(gl) - 1
    assert fig.axes[1].get_title(loc="left").startswith(f"{n:,} spacings")
    zero_marks = [ln for ln in fig.axes[0].lines if ln.get_marker() == "o"]
    assert sum(len(ln.get_xdata()) for ln in zero_marks) == int(
        np.sum((gz > 100) & (gz < 112)) + np.sum((gl > 100) & (gl < 112)))
