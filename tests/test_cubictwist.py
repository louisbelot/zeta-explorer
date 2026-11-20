import pytest

from zeta_explorer.cubic import OMEGA, Eis, lpolys, symbol
from zeta_explorer.cubictwist import (
    cubefree,
    family_size,
    weight,
)
from zeta_explorer.funcfield import Curve, legendre
from zeta_explorer.gf import deg, minimal_polys, monics, pdivmod, pgcd, pmul, primes_dividing

P7 = 7
IRR = Curve(P7, (), (), (1, 0, 1))
SPLIT = Curve(P7, (), (), (0, 1, 1))


def _brute_cubefree(p, N, avoid):
    out = []
    for n in range(N, 2 * N + 1):
        if n % 3:
            continue
        for D in monics(p, n):
            if deg(pgcd(D, avoid, p)):
                continue
            rad = primes_dividing(D, p)
            if sum(deg(R) for R in rad) != N:
                continue
            if any(not pdivmod(D, pmul(R, pmul(R, R, p), p), p)[1] for R in rad):
                continue
            out.append(D)
    return out


@pytest.mark.parametrize("N", [1, 2])
@pytest.mark.parametrize("E", [IRR, SPLIT])
def test_cubefree_family(N, E):
    assert sorted(cubefree(P7, N, E.disc)) == sorted(_brute_cubefree(P7, N, E.disc))


@pytest.mark.parametrize("N", [1, 2, 3, 4])
def test_family_size(N):
    for E in (IRR, SPLIT, legendre(P7)):
        assert family_size(E, N) == len(cubefree(P7, N, E.disc))


@pytest.mark.parametrize("N", [2, 3])
@pytest.mark.parametrize("d", [1, 2, 3])
def test_character_weight(N, d):
    E = SPLIT
    Ds = cubefree(P7, N, E.disc)
    bad = primes_dividing(E.disc, P7)
    Ps = [P for P in minimal_polys(P7, d)[1] if P not in bad][:6]
    for P, L in zip(Ps, lpolys(P7, Ps), strict=True):
        direct = Eis(0)
        for D in Ds:
            s = symbol(D, P, P7)
            if s is not None:
                direct = direct + OMEGA[s]
        badw = [(deg(R), symbol(R, P, P7)) for R in bad]
        assert direct == Eis(weight(P7, N, d, L, badw))
        coprime = sum(symbol(D, P, P7) is not None for D in Ds)
        assert coprime == weight(P7, N, d, None, [(deg(R), 0) for R in bad])
