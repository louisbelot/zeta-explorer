import mpmath
import numpy as np
import pytest

pytest.importorskip("flint")

from flint import acb  # noqa: E402

from zeta_explorer import skewes  # noqa: E402
from zeta_explorer.certified import precision  # noqa: E402
from zeta_explorer.mertens import arrays, certified_zeros  # noqa: E402

CHAO_PLYMEN = dict(omega=727.952074, eta=0.000198219, alpha=7.86e10, T=1131944.4718, A=7.79e6)


@pytest.fixture(scope="module")
def first_thousand():
    return arrays(certified_zeros(1420, bits=64))[0][:1000]


def test_chao_plymen_table_3(first_thousand):
    s = skewes.lehman_sum([437.7635, 437.782672], first_thousand)
    assert s == pytest.approx([0.88258736, 0.88873669], abs=2e-8)


def test_lehman_candidates(first_thousand):
    s = skewes.lehman_sum([727.952, 853.853, 2682.977], first_thousand)
    assert np.all(s > 0.95)
    assert np.all(skewes.lehman_sum([727.5, 853.5, 2682.5], first_thousand) < 0.5)


def test_gaussian_damps_the_sum(first_thousand):
    plain = skewes.lehman_sum(727.952, first_thousand)
    damped = skewes.lehman_sum(727.952, first_thousand, alpha=1e5)
    assert abs(damped[0]) < abs(plain[0])
    assert skewes.lehman_sum(727.952, first_thousand, alpha=1e12) == pytest.approx(plain, abs=1e-5)


def test_sum_is_vectorized_in_blocks(first_thousand, monkeypatch):
    w = np.linspace(700, 701, 7)
    full = skewes.lehman_sum(w, first_thousand)
    monkeypatch.setattr(skewes, "BLOCK", 10)
    assert skewes.lehman_sum(w, first_thousand) == pytest.approx(full, abs=1e-13)


def test_error_terms_match_chao_plymen():
    e = skewes.lehman_errors(**CHAO_PLYMEN)
    assert float(e["s2"]) == pytest.approx(5.9327e-50, rel=1e-4)
    assert float(e["s5"]) == pytest.approx(6.820541e-5, rel=1e-6)
    assert float(e["s6"]) == pytest.approx(9.75001e-5, rel=1e-5)
    assert e["s3"] < mpmath.mpf("1e-672") and e["s4"] < mpmath.mpf("1e-666")
    assert float(e["s1"]) == pytest.approx(2.1457 / 727.951875781, rel=1e-12)


def test_lehman_conditions():
    with pytest.raises(ValueError):
        skewes.lehman_errors(**{**CHAO_PLYMEN, "A": 7.8e6})
    with pytest.raises(ValueError):
        skewes.lehman_errors(**{**CHAO_PLYMEN, "T": 1e7})
    with pytest.raises(ValueError):
        skewes.lehman_errors(437.4, 0.000198219, 7.86e10, 1131944.4718, 7.79e6)
    assert skewes.lehman_errors(437.4, 0.000198219, 7.86e10, 1131944.4718, 7.79e6, False)


def test_zero_error_is_linear_in_delta(first_thousand):
    a = skewes.zero_error(727.952, 1e8, first_thousand, 1e-9)
    b = skewes.zero_error(727.952, 1e8, first_thousand, 2e-9)
    assert b == pytest.approx(2 * a, rel=1e-6)
    t = first_thousand[:50]
    h = 1e-7
    slope = (skewes.lehman_sum(727.952, t + h, 1e8) - skewes.lehman_sum(727.952, t, 1e8)) / h
    assert abs(slope[0]) * 1e-9 <= skewes.zero_error(727.952, 1e8, t, 1e-9) * 50


def test_crossing_interval_and_excess():
    c = skewes.Crossing(727.952074, 0.000198219, 7.86e10, 1131944.4718, 1.00319601829755,
                        2.05971416e-5, skewes.lehman_errors(**CHAO_PLYMEN))
    assert c.interval == ("1.39792136e+316", "1.39847567e+316")
    assert c.lower > 6e-5
    assert c.excess > 1e150


def test_polish_recovers_perturbed_zeros(first_thousand):
    t = first_thousand[600:]
    noisy = t + np.random.default_rng(0).uniform(-3e-8, 3e-8, len(t))
    assert np.max(np.abs(skewes.polish(noisy) - t)) < 1e-10
    with pytest.raises(ValueError):
        skewes.polish(t + 1e-3)


def test_bulk_zeros_keep_their_indices():
    g = skewes.bulk_zeros(45_000)
    assert len(g) == 45_000 and np.all(np.diff(g) > 0)
    with precision(128):
        for n in (1, 101, 102, 20_100, 20_101, 40_100, 45_000):
            assert abs(float(acb.zeta_zeros(n, 1)[0].imag.mid()) - g[n - 1]) < 1e-8
