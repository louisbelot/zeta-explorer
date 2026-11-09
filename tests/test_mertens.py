import math

import numpy as np
import pytest

pytest.importorskip("flint")

from flint import acb  # noqa: E402

from zeta_explorer import mertens as m  # noqa: E402
from zeta_explorer.certified import precision  # noqa: E402
from zeta_explorer.explicit import mobius  # noqa: E402

OEIS = [1, -1, 1, 2, -23, -48, 212, 1037, 1928, -222, -33722]


@pytest.fixture(scope="module")
def zeros():
    return m.certified_zeros(200, bits=256)


def test_sieve_matches_trial_division():
    mu = m.mobius_sieve(3000)
    assert mu[0] == 0
    assert all(mu[n] == mobius(n) for n in range(1, 3001))


def test_table_powers_of_ten():
    M = m.mertens_table(10**7)
    assert [int(M[10**k]) for k in range(8)] == OEIS[:8]


@pytest.mark.parametrize("k", [8, 9, 10])
def test_hyperbola_powers_of_ten(k):
    assert m.mertens(10**k) == OEIS[k]


def test_hyperbola_with_the_smallest_table():
    M = m.mertens_table(200_000)
    for x in [1, 2, 3, 17, 1000, 65536, 99991, 123457, 200_000]:
        assert m.mertens(x, u=math.isqrt(x) + 1) == M[x]


def test_zeros_are_complete_and_match_arb(zeros):
    assert len(zeros) == 79
    with precision(300):
        ref = acb.zeta_zeros(1, 79)
        assert all(z.gamma.overlaps(r.imag) for z, r in zip(zeros, ref, strict=True))
    assert max(float(z.gamma.rad()) for z in zeros) < 1e-70


def test_first_derivative_matches_mpmath(zeros):
    d = zeros[0].deriv
    assert abs(float(d.real.mid()) - 0.783296511867030928649657209239) < 1e-15
    assert abs(float(d.imag.mid()) - 0.124699829748171089409928491509) < 1e-15
    assert float(zeros[0].a.mid()) == pytest.approx(0.0891415213850345, rel=1e-12)


def test_kernel_is_nonnegative_and_vanishes_at_one():
    u = np.linspace(0, 1, 1001)
    f = m.kernel_np(u)
    assert f[0] == pytest.approx(1)
    assert abs(f[-1]) < 1e-15
    assert f.min() >= -1e-15


def test_explicit_formula_rebuilds_m(zeros):
    M = m.mertens_table(10**6)
    x = np.array([1000.5, 9999.5, 100_000.5, 523_456.5, 999_999.5])
    q = M[x.astype(int)] / np.sqrt(x)
    assert np.max(np.abs(m.q_from_zeros(x, zeros) - q)) < 0.04


def test_top_weights_match_odlyzko_te_riele():
    a = np.sort(m.arrays(m.certified_zeros(600, bits=64))[1][:400])[::-1]
    total = 2 * np.cumsum(a)
    assert total[24] == pytest.approx(0.780108, abs=5e-7)
    assert total[69] == pytest.approx(1.0787, abs=1e-4)
    assert np.argmax(total > 1) + 1 == 54
