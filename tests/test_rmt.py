import math

import numpy as np
import pytest

from zeta_explorer.rmt import (
    bessel,
    binned_chi2,
    char_poly_at_one,
    conditioned_n,
    even_orthogonal,
    even_orthogonal_n,
    excised_first,
    first_density,
    first_density_mix,
    fit_n,
    fredholm_det,
    gap,
    gap_chi2,
    gap_ks,
    independent_n,
    interaction_n,
    jacobi_angles,
    jacobi_moment,
    jacobi_n,
    odd_orthogonal_half,
    odd_orthogonal_n,
    rank_two_first,
    sine,
    small_value_law,
    so_angles,
    so_moment,
    symplectic,
    symplectic_n,
)

S = np.array([0.1, 0.4, 0.9])


def weyl_gap(alpha, density, nodes=80):
    u, w = np.polynomial.legendre.leggauss(nodes)

    def integral(lo):
        a = (np.pi - lo) / 2 * u + (np.pi + lo) / 2
        wa = w * (np.pi - lo) / 2
        t1, t2 = np.meshgrid(a, a)
        return np.sum(np.outer(wa, wa) * density(t1, t2))

    return integral(alpha) / integral(0.0)


def test_usp2_closed_form():
    k = symplectic_n(1)
    assert gap(k, S) == pytest.approx(1 - S + np.sin(2 * np.pi * S) / (2 * np.pi), abs=1e-14)
    assert first_density(k, S) == pytest.approx(1 - np.cos(2 * np.pi * S), abs=1e-14)


def test_so2_is_uniform():
    assert gap(even_orthogonal_n(1), S) == pytest.approx(1 - S, abs=1e-14)


@pytest.mark.parametrize("s", [0.3, 0.8, 1.5])
def test_rank_two_against_weyl_integral(s):
    alpha = np.pi * s / 2
    usp = weyl_gap(alpha, lambda a, b: (np.cos(a) - np.cos(b)) ** 2 * (np.sin(a) * np.sin(b)) ** 2)
    so = weyl_gap(alpha, lambda a, b: (np.cos(a) - np.cos(b)) ** 2)
    assert gap(symplectic_n(2), s) == pytest.approx(usp, abs=1e-13)
    assert gap(even_orthogonal_n(2), s) == pytest.approx(so, abs=1e-13)


def test_sine_gap_series():
    s = 0.1
    series = 1 - s + np.pi**2 * s**4 / 36 - np.pi**4 * s**6 / 675 + np.pi**6 * s**8 / 17640
    assert fredholm_det(sine, -s / 2, s / 2) == pytest.approx(series, abs=1e-11)


@pytest.mark.parametrize("s", [0.2, 0.7, 1.6])
def test_parity_split_of_the_sine_kernel(s):
    whole = fredholm_det(sine, -s, s)
    assert gap(symplectic, s) * gap(even_orthogonal, s) == pytest.approx(whole, abs=1e-14)


def test_first_zero_near_the_origin():
    assert first_density(symplectic, 1e-3) / 1e-6 == pytest.approx(2 * np.pi**2 / 3, rel=1e-5)
    assert first_density(even_orthogonal, 0.0) == 2


def test_first_density_is_a_density():
    s = np.linspace(0, 4, 2001)
    p = first_density(symplectic, s)
    assert np.all(p >= 0)
    assert np.trapezoid(p, s) == pytest.approx(1, abs=1e-6)
    assert gap(symplectic, 4.0) < 1e-12


def test_density_is_minus_gap_derivative():
    h = 1e-5
    for k in (symplectic, symplectic_n(3), even_orthogonal):
        fd = -(gap(k, 0.8 + h) - gap(k, 0.8 - h)) / (2 * h)
        assert first_density(k, 0.8) == pytest.approx(fd, abs=1e-8)


def test_finite_n_converges_like_one_over_n():
    limit = gap(symplectic, 1.0)
    err = [gap(symplectic_n(n), 1.0) - limit for n in (100, 1000)]
    assert err[0] / err[1] == pytest.approx(10, rel=0.02)


def test_quadrature_has_converged():
    assert gap(symplectic, 2.0, m=20) == pytest.approx(gap(symplectic, 2.0, m=60), abs=1e-14)


def test_fit_recovers_matrix_size():
    rng = np.random.default_rng(7)
    n = 3
    s = np.linspace(0, 1.6, 4001)
    cdf = 1 - gap(symplectic_n(n), s)
    sample = np.interp(rng.uniform(0, cdf[-1], 4000), cdf, s)
    best, lo, hi, _ = fit_n(sample, np.linspace(1.5, 8, 27))
    assert lo <= n <= hi
    assert hi - lo < 1.5
    assert math.isclose(best, n, abs_tol=0.6)


def test_kernels_reject_tiny_n():
    with pytest.raises(ValueError):
        symplectic_n(0.5)


def test_haar_so_one_level():
    n, size = 4, 20_000
    theta = so_angles(n, size, seed=3)
    counts, edges = np.histogram(theta.ravel(), bins=20, range=(0, np.pi))
    c = (edges[1:] + edges[:-1]) / 2
    width = edges[1] - edges[0]
    r1 = (2 * n - 1) / (2 * np.pi) + np.sin((2 * n - 1) * c) / (2 * np.pi * np.sin(c))
    z = (counts - size * width * r1) / np.sqrt(size * width * r1)
    assert np.abs(z).max() < 4.5
    assert theta.shape == (size, n)


def test_haar_so_first_eigenvalue_matches_fredholm():
    n = 3
    x = np.sort(n * so_angles(n, 20_000, seed=5)[:, 0] / np.pi)
    s = np.linspace(0, x[-1], 400)
    cdf = 1 - gap(even_orthogonal_n(n), s)
    ks = np.abs(np.interp(x, s, cdf) - np.arange(1, len(x) + 1) / len(x)).max()
    assert ks < 1.6 / np.sqrt(len(x))


def test_excised_ensemble_repels():
    plain = excised_first(4, 0.0, 10_000, seed=1)
    cut = excised_first(4, 0.05, 10_000, seed=1)
    assert len(plain) == 10_000
    assert 0 < len(cut) < len(plain)
    assert np.mean(cut) > np.mean(plain) + 0.05
    assert np.quantile(cut, 0.01) > 2 * np.quantile(plain, 0.01)


def test_fit_n_for_even_orthogonal():
    x = 3 * so_angles(3, 4000, seed=9)[:, 0] / np.pi
    best, lo, hi, _ = fit_n(x, np.linspace(2, 6, 17), family=even_orthogonal_n)
    assert lo <= 3 <= hi
    assert abs(best - 3) <= 0.5


def test_so3_closed_form():
    expect = 1 - S + np.sin(np.pi * S) / np.pi
    assert gap(odd_orthogonal_n(1), S) == pytest.approx(expect, abs=1e-14)


def test_odd_orthogonal_tends_to_symplectic_limit():
    assert gap(odd_orthogonal_n(40), S) == pytest.approx(gap(symplectic, S), abs=2e-3)


def test_haar_so_odd_one_level():
    n, size = 4, 20_000
    theta = so_angles(n, size, seed=4, odd=True)
    counts, edges = np.histogram(theta.ravel(), bins=20, range=(0, np.pi))
    u = np.linspace(0, np.pi, 20 * 64 + 1)[1:] - np.pi / (20 * 64 * 2)
    r1 = 2 * n / (2 * np.pi) - np.sin(2 * n * u) / (2 * np.pi * np.sin(u))
    expect = size * r1.reshape(20, 64).mean(axis=1) * (edges[1] - edges[0])
    z = (counts - expect) / np.sqrt(expect)
    assert np.abs(z).max() < 4.5
    assert theta.shape == (size, n)


def test_haar_so_odd_first_eigenvalue_matches_fredholm():
    n = 3
    x = np.sort(n * so_angles(n, 20_000, seed=6, odd=True)[:, 0] / np.pi)
    s = np.linspace(0, x[-1], 400)
    cdf = 1 - gap(odd_orthogonal_n(n), s)
    ks = np.abs(np.interp(x, s, cdf) - np.arange(1, len(x) + 1) / len(x)).max()
    assert ks < 1.6 / np.sqrt(len(x))


def test_so_moment_closed_forms():
    s = np.array([-0.4, 0.5, 1.0, 2.5])
    so3 = 4**s * np.array([math.gamma(1.5 + x) / (math.gamma(1.5) * math.gamma(2 + x)) for x in s])
    assert so_moment(1, s, odd=True) == pytest.approx(so3, rel=1e-13)
    assert so_moment(7, 1.0) == pytest.approx(2, rel=1e-13)
    assert so_moment(5, 0.0, odd=True) == 1


@pytest.mark.parametrize("odd", [False, True])
def test_so_moment_against_haar(odd):
    v = char_poly_at_one(so_angles(4, 40_000, seed=8, odd=odd))
    for s in (-0.3, 0.5, 1.0):
        mc = v**s
        assert abs(mc.mean() - so_moment(4, s, odd)) < 4 * mc.std() / np.sqrt(len(v))


def test_so_moment_rejects_the_pole():
    with pytest.raises(ValueError):
        so_moment(3, -1.5, odd=True)


def test_small_value_law():
    assert small_value_law(1, odd=True) == pytest.approx((1 / (4 * np.pi), 1.5), rel=1e-14)
    barnes_half = 0.603244281209446
    r, a = small_value_law(100)
    assert a == 0.5
    assert r == pytest.approx(2**-0.875 * barnes_half * np.pi**-0.25 * 100**0.375, rel=3e-3)


def test_small_value_law_against_haar():
    v = char_poly_at_one(so_angles(4, 60_000, seed=11, odd=True))
    r, a = small_value_law(4, odd=True)
    rho = 0.05
    k = int(np.sum(v < rho))
    assert abs(k - len(v) * r / a * rho**a) < 4 * np.sqrt(k)


def test_excised_odd_ensemble_repels():
    plain = excised_first(4, 0.0, 10_000, seed=2, odd=True)
    cut = excised_first(4, 2.0, 10_000, seed=2, odd=True)
    assert 0 < len(cut) < len(plain) == 10_000
    assert np.mean(cut) > np.mean(plain) + 0.03
    assert np.quantile(cut, 0.01) > 1.5 * np.quantile(plain, 0.01)


def test_binned_chi2():
    edges = np.linspace(0, 1, 5)
    x = np.array([0.1, 0.1, 0.6, 0.9])
    model = np.array([0.1, 0.3, 0.6, 0.9])
    assert binned_chi2(x, model, edges) == pytest.approx(2.0)
    assert binned_chi2(np.array([0.3]), np.array([0.1]), edges) == math.inf


def test_odd_orthogonal_half_unfolds_by_n_plus_half():
    n = 3
    x = np.sort((n + 0.5) * so_angles(n, 20_000, seed=12, odd=True)[:, 0] / np.pi)
    s = np.linspace(0, x[-1], 400)
    cdf = 1 - gap(odd_orthogonal_half(n), s)
    ks = np.abs(np.interp(x, s, cdf) - np.arange(1, len(x) + 1) / len(x)).max()
    assert ks < 1.6 / np.sqrt(len(x))


def test_odd_orthogonal_half_mean_density_is_one():
    k = odd_orthogonal_half(4)
    y = np.linspace(0, 4.5, 4001)
    assert np.trapezoid(k(y, y), y) == pytest.approx(4, abs=1e-6)


def test_first_density_mix():
    s = np.linspace(0, 1.5, 7)
    same = first_density_mix(np.full(10, 3.0), s, odd_orthogonal_half)
    assert same == pytest.approx(first_density(odd_orthogonal_half(3.0), s), abs=1e-14)
    two = first_density_mix([2.0, 2.0, 5.0], s, even_orthogonal_n, groups=2)
    expect = (2 * first_density(even_orthogonal_n(2.0), s)
              + first_density(even_orthogonal_n(5.0), s)) / 3
    assert two == pytest.approx(expect, abs=1e-14)


U = np.linspace(0.05, 0.9, 7)[:, None]


@pytest.mark.parametrize("n", [1, 3, 8])
def test_jacobi_reproduces_the_classical_kernels(n):
    for alpha, beta, k in ((-0.5, -0.5, even_orthogonal_n), (0.5, -0.5, odd_orthogonal_n),
                           (0.5, 0.5, symplectic_n)):
        assert jacobi_n(n, alpha, beta)(U, U.T) == pytest.approx(k(n)(U, U.T), abs=1e-13)


def test_independent_is_a_smaller_so_even():
    small = even_orthogonal_n(4)
    assert independent_n(5)(U, U.T) == pytest.approx(0.8 * small(0.8 * U, 0.8 * U.T), abs=1e-13)


def test_bessel_half_is_the_symplectic_limit():
    y = np.linspace(0.1, 2.7, 6)[None, :]
    assert bessel(0.5)(U, y) == pytest.approx(symplectic(U, y), abs=1e-14)
    assert bessel(0.5)(U, U) == pytest.approx(symplectic(U, U), abs=1e-14)


@pytest.mark.parametrize("alpha", [1.0, 1.5])
def test_jacobi_tends_to_bessel_like_one_over_n(alpha):
    limit = gap(bessel(alpha), 1.0)
    err = [gap(jacobi_n(n, alpha, scale=n), 1.0) - limit for n in (40, 400)]
    assert err[0] / err[1] == pytest.approx(10, rel=0.02)


def test_interaction_n2_closed_form():
    t = np.pi * S / 2
    inside = 1.5 * t - 2 * np.sin(t) + np.sin(2 * t) / 4
    assert gap(interaction_n(2), S) == pytest.approx(1 - inside / (1.5 * np.pi), abs=1e-14)


def test_conditioned_n2_closed_form():
    def primitive(phi):
        return -np.cos(phi) + np.cos(phi) ** 3 / 3

    t = np.pi * S / 4
    expect = (primitive(np.pi / 2) - primitive(t)) / (primitive(np.pi / 2) - primitive(0))
    assert gap(conditioned_n(2), S) == pytest.approx(expect, abs=1e-13)


@pytest.mark.parametrize("s", [0.4, 1.0, 1.7])
def test_interaction_n3_against_weyl_integral(s):
    weyl = weyl_gap(np.pi * s / 3, lambda a, b: (np.cos(a) - np.cos(b)) ** 2
                    * ((1 - np.cos(a)) * (1 - np.cos(b))) ** 2)
    assert gap(interaction_n(3), s) == pytest.approx(weyl, abs=1e-13)


def test_forced_pair_repels_like_a_power():
    for alpha, k in ((1.5, interaction_n(8)), (1.0, conditioned_n(8))):
        ratio = first_density(k, 0.02) / first_density(k, 0.01)
        assert np.log2(ratio) == pytest.approx(2 * alpha + 1, abs=0.01)


@pytest.mark.parametrize("model", ["independent", "interaction", "conditioned"])
def test_rank_two_monte_carlo_matches_fredholm(model):
    kernel = {"independent": independent_n, "interaction": interaction_n,
              "conditioned": conditioned_n}[model](4)
    x = rank_two_first(4, model, 20_000, seed=5)
    assert gap_ks(x, kernel) < 1.6 / np.sqrt(len(x))


def test_jacobi_moment_closed_form():
    s = np.array([-1.2, 0.5, 2.0])
    expect = 4**s * np.array([math.gamma(2.5 + x) * 2 / (math.gamma(2.5) * math.gamma(3 + x))
                              for x in s])
    assert jacobi_moment(1, s, 1.5) == pytest.approx(expect, rel=1e-13)
    assert jacobi_moment(4, 0.5, 0.5) == so_moment(4, 0.5, odd=True)


def test_jacobi_angles_reweighting():
    v = char_poly_at_one(jacobi_angles(3, 1.5, 40_000, seed=2))
    for s in (0.5, 1.0):
        mc = v**s
        assert abs(mc.mean() - jacobi_moment(3, s, 1.5)) < 4 * mc.std() / np.sqrt(len(v))
    with pytest.raises(ValueError):
        jacobi_angles(3, 0.2, 10)


def test_rank_two_cutoff_removes_small_values():
    plain = rank_two_first(4, "interaction", 10_000, seed=4)
    cut = rank_two_first(4, "interaction", 10_000, seed=4, cutoff=1.0)
    assert 0 < len(cut) < len(plain) == 10_000
    assert np.mean(cut) > np.mean(plain)


def test_gap_chi2_and_ks():
    flat = even_orthogonal_n(1)
    edges = np.array([0.0, 0.5, 1.0])
    assert gap_chi2(np.linspace(0.05, 0.95, 10), flat, edges) == pytest.approx(0, abs=1e-20)
    assert gap_chi2(np.full(10, 0.2), flat, edges) == pytest.approx(10)
    assert gap_ks(np.array([0.5, 1.0]), flat) == pytest.approx(0.5)
