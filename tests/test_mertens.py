import math

import pytest

pytest.importorskip("flint")

from zeta_explorer import mertens as m  # noqa: E402
from zeta_explorer.explicit import mobius  # noqa: E402

OEIS = [1, -1, 1, 2, -23, -48, 212, 1037, 1928, -222, -33722]


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
