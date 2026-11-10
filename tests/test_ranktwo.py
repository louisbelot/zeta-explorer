import math

import numpy as np
import pytest

from zeta_explorer.elliptic import Curve, EllipticL, low_zeros
from zeta_explorer.ranktwo import compare, family, first_zeros, matrix_size, rank_two
from zeta_explorer.rmt import rank_two_first

E11 = Curve(0, -1, 1, 0, 0)


def test_rank_two_twist():
    L = EllipticL(E11, 232)
    t = rank_two(L)
    assert t.d == 232
    assert t.log_A == pytest.approx(math.log(L.A))
    assert t.second > 0
    assert t.x == pytest.approx(low_zeros(L).x)
    assert 0 < t.x[0] < t.x[1]


def test_rank_two_rejects_rank_zero():
    with pytest.raises(ValueError):
        rank_two(EllipticL(E11, 5))


def test_family_finds_the_waldspurger_zeros():
    fam = family(E11, 5, 400, workers=1)
    assert [t.d for t in fam] == [232, 265, 273, 344, 364]
    assert matrix_size(fam) == pytest.approx(np.mean([math.log(EllipticL(E11, t.d).A)
                                                      for t in fam]))
    assert first_zeros(fam) == pytest.approx([t.x[0] for t in fam])


def test_compare_picks_the_sampled_model():
    x = rank_two_first(5, "interaction", 3000, seed=7)
    table = compare(x, 5, np.linspace(0, 2.5, 11), grid=range(3, 8))
    ks = {name: row[1] for name, row in table.items()}
    assert min(ks, key=ks.get) == "interaction"
    assert table["interaction"][2] == 5
