import itertools

import numpy as np
import pytest

from zeta_explorer.cubic import (
    OMEGA,
    Eis,
    index,
    symbol,
    zeta3,
)
from zeta_explorer.gf import field, minimal_polys, pmul

P7 = 7


def test_eisenstein_ring():
    w = OMEGA[1]
    assert w**3 == Eis(1)
    assert Eis(1) + w + w * w == Eis(0)
    x, y = Eis(3, -2), Eis(-1, 5)
    assert (x * y).norm == x.norm * y.norm
    assert (x * y).conj() == x.conj() * y.conj()
    assert complex(x * y) == pytest.approx(complex(x) * complex(y))
    assert (x * x.conj()) == Eis(x.norm)
    assert x.tr == pytest.approx(2 * complex(x).real)


def test_zeta3():
    assert zeta3(7) == 2 and zeta3(13) == 3
    with pytest.raises(ValueError):
        zeta3(5)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_symbol_matches_field_index(d):
    F = field(P7, d)
    alpha, Ps = minimal_polys(P7, d)
    idx = index(F)
    rng = np.random.default_rng(d)
    for i in rng.choice(len(Ps), size=min(5, len(Ps)), replace=False):
        for A in [tuple(int(c) for c in rng.integers(0, P7, size=4)) for _ in range(6)]:
            v = idx[F.polyval(A, alpha[i:i + 1])][0]
            assert symbol(A, Ps[i], P7) == (None if v < 0 else v)


def test_symbol_is_multiplicative():
    P = minimal_polys(P7, 2)[1][3]
    for A, B in itertools.product([(1, 1), (3, 0, 1), (2, 5)], [(4, 1), (1, 2, 3)]):
        assert symbol(pmul(A, B, P7), P, P7) == (symbol(A, P, P7) + symbol(B, P, P7)) % 3


@pytest.mark.parametrize("da,db", [(1, 2), (2, 2), (2, 3), (1, 3)])
def test_cubic_reciprocity(da, db):
    As = minimal_polys(P7, da)[1][:12]
    Bs = minimal_polys(P7, db)[1][:12]
    for A, B in itertools.product(As, Bs):
        if A != B:
            assert symbol(A, B, P7) == symbol(B, A, P7)
