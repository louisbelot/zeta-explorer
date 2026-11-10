import math

import numpy as np
import pytest

from zeta_explorer.elliptic import Curve, EllipticL, low_zeros
from zeta_explorer.ranktwo import (
    compare,
    excision_scan,
    family,
    first_zeros,
    ks2,
    matrix_size,
    rank_two,
)
from zeta_explorer.rmt import char_poly_at_one, jacobi_angles, rank_two_first

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


def test_ks2():
    assert ks2(np.array([1.0, 2.0]), np.array([1.0, 2.0])) == 0
    assert ks2(np.array([0.0, 1.0]), np.array([2.0, 3.0, 4.0])) == 1
    assert ks2(np.array([1.0, 3.0]), np.array([2.0, 4.0])) == 0.5


def test_excision_scan_finds_the_planted_cut():
    theta = jacobi_angles(4, 1.5, 40_000, seed=11)
    v = char_poly_at_one(theta)
    x = 5 * theta[v >= np.quantile(v, 0.4), 0][:4000] / np.pi
    fractions = np.linspace(0, 0.8, 9)
    cuts, ks = excision_scan(x, 5, fractions, size=40_000, seed=12)
    assert np.all(np.diff(cuts) > 0)
    assert fractions[np.argmin(ks)] == pytest.approx(0.4)
    assert ks.min() < 1.36 / np.sqrt(4000)
