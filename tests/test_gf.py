import math

import numpy as np
import pytest

from zeta_explorer.gf import (
    GF,
    correlate,
    cubic_sums,
    distinct_roots,
    factor,
    field,
    monic,
    monics,
    padd,
    pdivmod,
    peval,
    pgcd,
    pmul,
    ppowmod,
    primes_dividing,
    squarefree,
)


@pytest.fixture(params=[(3, 4), (5, 3), (7, 2), (11, 1)])
def F(request):
    return field(*request.param)


def test_field_axioms(F):
    rng = np.random.default_rng(1)
    a, b, c = rng.integers(0, F.q, (3, 500))
    assert np.array_equal(F.mul(a, F.add(b, c)), F.add(F.mul(a, b), F.mul(a, c)))
    assert np.array_equal(F.add(a, F.neg(a)), np.zeros_like(a))
    nz = a[a > 0]
    assert np.all(F.mul(nz, F.power(nz, F.q - 2)) == 1)


def test_quadratic_character(F):
    x = np.arange(1, F.q)
    assert np.sum(F.chi == 1) == (F.q - 1) // 2
    assert np.array_equal(F.chi[F.mul(x, x[::-1])], F.chi[x] * F.chi[x[::-1]])
    assert F.chi[F.neg(1)] == (-1) ** ((F.q - 1) // 2)


def test_frobenius(F):
    x = np.arange(F.q)
    y = x[::-1]
    assert np.array_equal(F.frobenius(F.add(x, y)), F.add(F.frobenius(x), F.frobenius(y)))
    assert np.array_equal(np.flatnonzero(F.frobenius(x) == x), np.arange(F.p))
    cur = x
    for _ in range(F.m):
        cur = F.frobenius(cur)
    assert np.array_equal(cur, x)


def test_orbit_count(F):
    reps, inv = F.orbits
    necklaces = sum(F.p ** math.gcd(F.m, k) for k in range(F.m)) // F.m
    assert len(reps) == necklaces
    assert np.array_equal(reps[inv][F.frobenius(np.arange(F.q))], reps[inv])


def test_polyval_on_the_prime_field():
    F = field(13, 1)
    c = (3, 0, 5, 1, 12)
    assert F.polyval(c, np.arange(13)).tolist() == [peval(c, x, 13) for x in range(13)]


def test_cubic_sums_against_legendre_symbols():
    p = 11
    F = field(p, 1)
    chi = [0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1 for x in range(1, p)]
    for A, B, C in [(0, 1, 3), (4, 0, 7), (2, 9, 0)]:
        brute = sum(chi[(x**3 + A * x * x + B * x + C) % p] for x in range(p))
        assert cubic_sums(F, A, B, C)[0] == brute


def test_cubic_sums_follow_frobenius_powers():
    p = 7
    for A, B, C in [(0, 1, 3), (1, 0, 2), (3, 5, 6)]:
        a1 = -cubic_sums(field(p, 1), A, B, C)[0]
        a2 = -cubic_sums(field(p, 2), A, B, C)[0]
        a3 = -cubic_sums(field(p, 3), A, B, C)[0]
        assert a2 == a1 * a1 - 2 * p
        assert a3 == a1**3 - 3 * p * a1


def test_correlate_matches_the_direct_sum():
    F = field(3, 3)
    rng = np.random.default_rng(2)
    g, h = rng.integers(-3, 4, (2, F.q))
    x = np.arange(F.q)
    direct = [int(np.sum(g * h[F.add(x, F.neg(np.full_like(x, t)))])) for t in range(F.q)]
    assert correlate(F, g, h).tolist() == direct


def test_polynomial_division_and_gcd():
    p = 7
    a, b = (3, 1, 4, 1, 5, 2), (2, 6, 1)
    q, r = pdivmod(a, b, p)
    assert len(r) < len(b)
    assert padd(pmul(q, b, p), r, p) == a
    c = (1, 1)
    assert pgcd(pmul(a, c, p), pmul(b, c, p), p) == monic(pmul(pgcd(a, b, p), c, p), p)
    assert ppowmod((0, 1), p, (1, 0, 0, 1), p) == pdivmod((0,) * p + (1,), (1, 0, 0, 1), p)[1]


def test_distinct_roots():
    p = 5
    cube = pmul(pmul((-1, 1), (-1, 1), p), (-1, 1), p)
    assert distinct_roots(pmul(cube, (-2, 1), p), p) == 2
    assert distinct_roots((0, -1, 0, 0, 0, 1), p) == 5
    irr = (-2, 0, 1)
    assert distinct_roots(pmul(irr, irr, p), p) == 2
    six = (1,)
    for _ in range(6):
        six = pmul(six, (1, 0, 1), 3)
    assert distinct_roots(pmul(six, (0, 1), 3), 3) == 3
    assert not squarefree(pmul(irr, irr, p), p)
    assert squarefree(pmul(irr, (0, 1), p), p)


def _mobius(n):
    out, d = 1, 2
    while d * d <= n:
        if n % d == 0:
            n //= d
            if n % d == 0:
                return 0
            out = -out
        d += 1
    return -out if n > 1 else out


@pytest.mark.parametrize("p,n", [(3, 1), (3, 2), (3, 4), (5, 3), (7, 2)])
def test_factor_counts_irreducibles(p, n):
    irreducible = 0
    for f in monics(p, n):
        if not squarefree(f, p):
            continue
        parts = factor(f, p)
        prod = (1,)
        for P in parts:
            prod = pmul(prod, P, p)
        assert prod == f
        irreducible += len(parts) == 1
    assert irreducible == sum(_mobius(n // d) * p**d for d in range(1, n + 1) if n % d == 0) // n


def test_factor_rejects_repeated_factors():
    with pytest.raises(ValueError):
        factor(pmul((1, 1), (1, 1), 5), 5)


def test_primes_dividing_survives_pth_powers():
    p, P = 3, (1, 0, 1)
    cube = pmul(pmul(P, P, p), P, p)
    assert primes_dividing(cube, p) == [P]
    f = pmul(pmul(cube, (0, 1), p), pmul((1, 1), (1, 1), p), p)
    assert primes_dividing(f, p) == [(0, 1), (1, 1), P]


def test_rejects_even_and_composite():
    with pytest.raises(ValueError):
        GF(2, 3)
    with pytest.raises(ValueError):
        GF(9, 1)


def test_degrees_count_elements_of_each_degree(F):
    size = F.degrees(np.arange(F.q))
    for d in range(1, F.m + 1):
        want = sum(_mobius(d // e) * F.p**e for e in range(1, d + 1) if d % e == 0)
        assert np.sum(size == d) == (want if F.m % d == 0 else 0)
