import numpy as np
import pytest

from zeta_explorer.cubic import OMEGA, Eis, logderiv, lpolys, symbol
from zeta_explorer.cubictwist import (
    cubefree,
    cubic_sums,
    cubic_twist,
    exact_traces,
    family_size,
    hecke_split,
    j0_frobenius,
    mean_traces,
    orthogonal_family,
    unitary_family,
    unitary_sums,
    weight,
)
from zeta_explorer.funcfield import Curve, legendre, lpoly
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


@pytest.mark.parametrize("E", [IRR, SPLIT])
def test_orthogonal_family_from_point_counts(E):
    fam = orthogonal_family(E, 2)
    for t in fam[::7]:
        assert lpoly(cubic_twist(E, t.D)) == t.L
    for t in fam:
        assert t.L.N == E.degree + 4
        r = t.L.roots()
        assert np.allclose(np.abs(r), 1)


def test_unitary_family_is_a_polynomial_on_the_circle():
    E = legendre(P7)
    fam = unitary_family(E, 2)
    M = E.degree + 4
    extra = unitary_sums(E, [t.D for t in fam], M + 1)
    for t, s in zip(fam, extra, strict=True):
        assert np.allclose(np.abs(t.roots()), 1)
        assert t.c[M].norm == P7 ** (2 * M)
        assert logderiv(list(t.c), M + 1)[M] == s


@pytest.mark.parametrize("N", [2, 3])
@pytest.mark.parametrize("E", [IRR, SPLIT])
def test_exact_traces_orthogonal(N, E):
    Ds = cubefree(P7, N, E.disc)
    brute = mean_traces([list(cubic_sums(E, Ds, m)) for m in range(1, 5)], P7)
    assert exact_traces(E, N, 4) == brute


@pytest.mark.parametrize("N", [2, 3])
def test_exact_traces_unitary(N):
    E = legendre(P7)
    Ds = cubefree(P7, N, E.disc)
    brute = mean_traces([unitary_sums(E, Ds, m) for m in range(1, 5)], P7)
    assert exact_traces(E, N, 4) == brute


def test_both_twists_of_a_j0_curve_share_means():
    Ds = cubefree(P7, 2, SPLIT.disc)
    a = mean_traces([list(cubic_sums(SPLIT, Ds, m)) for m in range(1, 5)], P7)
    b = mean_traces([unitary_sums(SPLIT, Ds, m) for m in range(1, 5)], P7)
    assert a == b


def test_scaling_kills_means_for_linear_B():
    E = Curve(P7, (), (), (0, 1))
    assert all(x == 0 for x in exact_traces(E, 3, 4))
    Ds = cubefree(P7, 2, E.disc)
    assert all(x == 0 for x in mean_traces([list(cubic_sums(E, Ds, m)) for m in (1, 2)], P7))


def test_diagonal_is_most_of_the_mean_at_large_N():
    full = exact_traces(SPLIT, 8, 3)
    diag = exact_traces(SPLIT, 8, 3, offdiagonal=False)
    assert diag[0] == diag[1] == 0
    assert all(abs(float(f - g)) < 2e-3 for f, g in zip(full, diag, strict=True))


def test_cubic_twist_needs_j0():
    with pytest.raises(ValueError):
        cubic_twist(legendre(P7), (1, 1))


@pytest.mark.parametrize("p", [7, 13, 19])
def test_j0_frobenius_has_norm_p(p):
    assert j0_frobenius(p).norm == p


@pytest.mark.parametrize("E", [IRR, SPLIT])
def test_twists_split_into_conjugate_hecke_l_functions(E):
    pi = j0_frobenius(P7)
    for t in orthogonal_family(E, 2):
        H = hecke_split(E, t.D)
        a = [h * pi**k for k, h in enumerate(H)]
        b = [x.conj() for x in a]
        h = len(H) - 1
        prod = [sum((a[i] * b[k - i] for i in range(max(0, k - h), min(k, h) + 1)), Eis(0))
                for k in range(2 * h + 1)]
        assert prod == [Eis(c) for c in t.L.c]
