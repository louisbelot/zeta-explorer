import math

import pytest

from zeta_explorer.siegel import (
    bound,
    constants,
)

TABLE = {1.6: (0.22675, 1.1614), 1.3: (0.23083, 1.1805), 1.2: (0.23222, 1.1870),
         1.1: (0.23362, 1.1936)}


def r_of(lam):
    return lam / (10 * math.log(10))


@pytest.mark.parametrize("lam", sorted(TABLE))
def test_constants_match_table(lam):
    k = constants(r_of(lam))
    phi, E = TABLE[lam]
    assert math.ceil(k.phi * 1e5) / 1e5 == pytest.approx(phi, abs=1e-12)
    assert math.ceil(k.E * 1e4) / 1e4 == pytest.approx(E, abs=1e-12)


def test_phi_tends_to_quarter():
    assert constants(1e-6).phi == pytest.approx(0.25, abs=1e-5)
    rs = [0.02, 0.05, 0.1, 0.2]
    phis = [constants(r).phi for r in rs]
    assert phis == sorted(phis, reverse=True)


def test_bound_is_theorem_at_the_edge():
    k = constants(r_of(1.3))
    R = k.r + 7 / 8
    for q in (10**5, 10**8, 10**10):
        L = math.log(q)
        beta = 1 - 0.2 / L
        gap = 1 + k.r - beta
        want = 1 / k.r - 1 / gap + gap / R**2 + k.phi * L + k.E
        assert float(bound(q, k)) == pytest.approx(want, rel=1e-13)
