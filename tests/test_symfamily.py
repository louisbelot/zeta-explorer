import math

import mpmath
import numpy as np
import pytest

from zeta_explorer.hecke import dimension, eigenforms, petersson, trace
from zeta_explorer.mellin import MellinLine
from zeta_explorer.selfdual import turing
from zeta_explorer.symfamily import (
    Family,
    count_zeros,
    family,
    harmonic_weight,
    lambda_power,
    low_zeros,
    plancherel_mean,
    sym_square,
    unfold,
)
from zeta_explorer.sympower import delta_power
from zeta_explorer.sympower import petersson as delta_norm


def test_weight_12_is_sym_square_of_delta():
    L = sym_square(eigenforms(12)[0])
    ref = delta_power(2)
    for s in (0.5 + 3j, 0.8 + 11j, 1.0):
        assert L(s) == pytest.approx(ref(s), rel=1e-12)


def test_harmonic_weight_of_delta_from_its_petersson_norm():
    k = 12
    w = mpmath.gamma(k - 1) / ((4 * mpmath.pi) ** (k - 1) * delta_norm())
    assert harmonic_weight(eigenforms(k)[0]) == pytest.approx(float(w), rel=1e-12)


@pytest.mark.parametrize("k", [24, 30, 36, 44])
def test_petersson_trace_formula(k):
    F = eigenforms(k)
    w = np.array([harmonic_weight(f) for f in F])
    lam = {1: np.ones(len(F))}
    for n, i in ((2, 0), (3, 1), (5, 2)):
        lam[n] = np.array([f.lam[i] for f in F])
    lam[4] = lambda_power(2, lam[2])
    for m, n in [(1, 1), (1, 2), (2, 2), (2, 3), (3, 3), (1, 4), (4, 5), (5, 5)]:
        assert np.sum(w * lam[m] * lam[n]) == pytest.approx(petersson(k, m, n), abs=1e-11)


def test_lambda_power_is_the_hecke_recursion():
    f = eigenforms(24)[1]
    x = f.lam[0]
    assert lambda_power(0, x) == 1 and lambda_power(1, x) == x
    theta = f.angles()[0]
    for n in range(2, 8):
        assert lambda_power(n, x) == pytest.approx(math.sin((n + 1) * theta) / math.sin(theta))


def _natural_mean(p, K):
    ks = range(K, 2 * K, 2)
    return sum(trace(k, p * p) / p ** (k - 1) for k in ks) / sum(dimension(k) for k in ks)


def test_natural_mean_of_lambda_p2_tends_to_plancherel():
    for p in (2, 3, 5, 7):
        gap = [_natural_mean(p, K) - plancherel_mean(p, 2) for K in (100, 400)]
        assert -0.04 < gap[0] < 0 and 0.15 < gap[1] / gap[0] < 0.35
    assert plancherel_mean(3, 3) == 0 and plancherel_mean(3, 4) == 1 / 9


def test_harmonic_mean_of_lambda_p2_vanishes_below_the_bessel_transition():
    F = eigenforms(60)
    w = np.array([harmonic_weight(f) for f in F])
    for i, p in enumerate((2, 3, 5, 7)):
        mean = np.sum(w * lambda_power(2, np.array([f.lam[i] for f in F])))
        assert mean == pytest.approx(petersson(60, 1, p * p), abs=1e-11)
        if 4 * math.pi * p < 59:
            assert abs(mean) < 1e-6
        else:
            assert abs(mean) > 0.1


def test_zeros_agree_with_the_pointwise_afe_and_turing():
    f = eigenforms(36)[2]
    L = sym_square(f)
    x = low_zeros(L, 6.0)
    T = 1.0
    while unfold(L, T) < 6.5:
        T *= 1.5
    direct = L.zeros(T)
    direct = direct[unfold(L, direct) <= 6.0]
    assert len(x) == len(direct) >= 4
    assert unfold(L, direct) == pytest.approx(x, abs=1e-9)
    assert np.all(np.diff(x) > 0) and x[0] > 0
    full = MellinLine(L, 60.0).zeros()
    assert turing(L, full, 42.0, 60.0).complete


def test_family_density_and_weights():
    F = family(range(12, 52, 2), x_max=3.0)
    assert F.size == sum(dimension(k) for k in range(12, 52, 2))
    assert np.all(F.weight > 0)
    big = F.select(F.k >= 30)
    for k in set(big.k.tolist()):
        assert big.weight[big.k == k].sum() == pytest.approx(1.0, abs=1e-3)
    x, rho = F.density(bins=10)
    assert len(x) == 10 and np.all(rho >= 0)
    count = sum(len(z) for z in F.zeros)
    assert rho.sum() * 0.3 * F.size == pytest.approx(count)
    _, rho_h = F.density(bins=10, harmonic=True)
    assert rho_h.sum() * 0.3 == pytest.approx(
        np.sum(F.weight * [len(z) for z in F.zeros]) / F.weight.sum())
    assert np.all(F.first() > 0) and np.all(F.first() <= 3.0)


def test_chi2_by_hand_and_round_trip(tmp_path):
    F = Family(np.array([12, 12]), np.array([1.0, 3.0]),
               [np.array([0.2, 0.7]), np.array([0.3])], 1.0)
    flat = lambda x: np.ones_like(x)  # noqa: E731
    assert F.chi2(flat, bins=2, x_max=1.0) == pytest.approx(1.0)
    assert F.chi2(lambda x: 2 * x, bins=2, x_max=1.0) == pytest.approx(1.5**2 / 0.5 + 0.5**2 / 1.5)
    assert F.chi2(flat, bins=2, x_max=1.0, harmonic=True) == pytest.approx(1.25 / 1.1)
    F.save(tmp_path / "fam.npz")
    G = Family.load(tmp_path / "fam.npz")
    assert G.k.tolist() == F.k.tolist() and G.weight.tolist() == F.weight.tolist()
    assert [z.tolist() for z in G.zeros] == [z.tolist() for z in F.zeros] and G.x_max == 1.0


def test_argument_principle_count_matches_the_scan():
    for k, i in ((12, 0), (60, 3), (120, 7)):
        L = sym_square(eigenforms(k)[i])
        z = MellinLine(L, 40.0).zeros(per_gap=32)
        for j in (3, 6):
            assert count_zeros(L, (z[j - 1] + z[j]) / 2) == j
    with pytest.raises(ArithmeticError):
        count_zeros(L, z[2] + 1e-9)
