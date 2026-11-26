import math
from fractions import Fraction

import numpy as np
import pytest

from zeta_explorer.brauer import ratio, scan
from zeta_explorer.bsd import arith
from zeta_explorer.funcfield import family, legendre
from zeta_explorer.heights import rank_one


@pytest.fixture(scope="module")
def fam5():
    return family(legendre(5), 2)


def test_rank_zero_is_sha(fam5):
    for t in fam5:
        if t.rank == 0:
            r = ratio(t)
            assert r.sha_reg == arith(t).sha
            assert r.bs == pytest.approx(math.log(arith(t).sha) / (r.chi * math.log(5)))


def test_rank_one_is_sha_times_height(fam5):
    for t in fam5:
        if t.rank == 1:
            r, s = ratio(t), rank_one(t)
            assert r.sha_reg == s.sha * s.generator.height


def test_rational_surface():
    for t in family(legendre(7), 1):
        r = ratio(t)
        assert r.chi == 1
        assert r.sha_reg == (1 if t.rank == 0 else Fraction(1, 2))


def test_ratio_from_parts(fam5):
    for t in fam5:
        r = ratio(t)
        expect = (math.log(r.lstar) + (r.chi - 1) * math.log(5) + math.log(16)
                  - math.log(r.tamagawa)) / (r.chi * math.log(5))
        assert r.bs == pytest.approx(expect)


def test_scan():
    S = scan(3, 4)
    assert S.chi == 3 and len(S.bs) == 29
    assert np.all(np.isin(S.rank, [0, 1, 2]))
    assert S.bs == pytest.approx(np.log(S.sha_reg) / S.log_h)
    assert 0 < S.trivial < 1
