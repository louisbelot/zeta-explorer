import mpmath
import pytest

from zeta_explorer import barrier as P

N0, N1 = 69_098, 1_500_000


def test_erfi_tail_is_the_integral():
    sigma, t = 1.3, 0.2
    exact = mpmath.quad(lambda a: mpmath.exp(t / 4 * mpmath.log(a) ** 2) / a**sigma, [100, 5000])
    assert P.erfi_tail(100, 5000, sigma, t) == pytest.approx(float(exact), rel=1e-12)


def test_tails_bound_the_sum():
    sigma, t, M, N = 1.3, 0.2, 100, 5000
    total = P.bsum(N, sigma, t)
    head = P.bsum(M, sigma, t)
    assert total <= head + P.erfi_tail(M, N, sigma, t) <= head + P.lemma_tail(M, N, sigma, t)
    assert total - head > P.erfi_tail(M + 1, N, sigma, t)


def test_lemma_needs_room():
    with pytest.raises(ValueError):
        P.lemma_tail(10, 1e6, 0.5, 0.2)


def test_far_region_from_the_paper():
    A, B = P.far_bound(N1, 0.2, 0.2, N0)
    assert B == pytest.approx(0.0740, abs=1e-4)
    assert A == pytest.approx(1.906, abs=1e-3)
    assert P.far_sup(N1, 1e30, 0.2, 0.2, N0) == pytest.approx(A + B, rel=1e-9)
    assert A + B < 2 - 1.25e-3


def test_top_edge_from_the_paper():
    assert 0.6 < P.top_bound(N0, N1, 0.2) < 0.8
