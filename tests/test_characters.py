import math

import numpy as np
import pytest

from zeta_explorer.characters import (
    characters,
    class_number,
    conrey,
    conrey_generator,
    factor,
    fundamental_unit,
    is_fundamental,
    kronecker,
    primitive_characters,
    quadratic,
    totient,
)

MODULI = [1, 2, 3, 4, 5, 7, 8, 9, 12, 15, 16, 20, 24, 25, 27, 32, 45, 63, 64, 100]


def mobius(n):
    f = factor(n)
    return 0 if any(e > 1 for _, e in f) else (-1) ** len(f)


def test_factor():
    assert factor(1) == []
    assert factor(360) == [(2, 3), (3, 2), (5, 1)]
    assert factor(9973) == [(9973, 1)]


@pytest.mark.parametrize("p, g", [(3, 2), (5, 2), (7, 3), (11, 2), (13, 2), (23, 5), (41, 6)])
def test_conrey_generator_is_least_primitive_root(p, g):
    assert conrey_generator(p) == g


def test_conrey_generator_skips_roots_that_fail_mod_p_squared():
    assert conrey_generator(40487) == 10


@pytest.mark.parametrize("q", MODULI)
def test_orthogonality(q):
    chars = characters(q)
    assert len(chars) == totient(q)
    table = np.array([c.values for c in chars])
    gram = table @ table.conj().T
    assert np.allclose(gram, totient(q) * np.eye(len(chars)), atol=1e-9)
    coprime = np.array([math.gcd(m, q) == 1 for m in range(q)])
    cols = table.conj().T @ table
    assert np.allclose(cols, totient(q) * np.diag(coprime), atol=1e-9)


@pytest.mark.parametrize("q", MODULI)
def test_number_of_primitive_characters(q):
    expected = sum(mobius(q // d) * totient(d) for d in range(1, q + 1) if q % d == 0)
    assert len(primitive_characters(q)) == expected


@pytest.mark.parametrize("q", [7, 16, 24, 45, 63])
def test_conrey_labels_multiply_and_are_symmetric(q):
    units = [n for n in range(1, q) if math.gcd(n, q) == 1]
    for m in units:
        for n in units:
            assert np.allclose((conrey(q, m) * conrey(q, n)).values, conrey(q, m * n % q).values)
            assert conrey(q, m)(n) == pytest.approx(conrey(q, n)(m), abs=1e-12)


@pytest.mark.parametrize("q", [9, 12, 16, 20, 27, 40, 63, 100])
def test_primitive_induces_the_character(q):
    for chi in characters(q):
        prim = chi.primitive()
        assert prim.is_primitive and prim.q == chi.conductor
        for m in range(q):
            if math.gcd(m, q) == 1:
                assert chi(m) == pytest.approx(prim(m), abs=1e-12)


@pytest.mark.parametrize(
    "q, n, order, parity, conductor",
    [(7, 6, 2, 1, 7), (7, 2, 3, 0, 7), (7, 3, 6, 1, 7), (8, 3, 2, 1, 8), (8, 5, 2, 0, 8),
     (5, 2, 4, 1, 5), (9, 4, 3, 0, 9), (12, 5, 2, 1, 3), (16, 9, 2, 0, 8), (32, 17, 2, 0, 8)],
)
def test_lmfdb_character_data(q, n, order, parity, conductor):
    chi = conrey(q, n)
    assert (chi.order, chi.parity, chi.conductor) == (order, parity, conductor)


def test_conrey_values_mod_seven():
    chi = conrey(7, 3)
    assert chi(3) == pytest.approx(complex(0.5, math.sqrt(3) / 2))
    assert chi(6) == -1
    assert chi(14) == 0


@pytest.mark.parametrize("q", [3, 4, 5, 8, 11, 12, 13, 15, 16, 21, 25, 27, 40, 99])
def test_gauss_sums_of_primitive_characters(q):
    for _, chi in primitive_characters(q):
        assert abs(chi.gauss_sum()) == pytest.approx(math.sqrt(q))
        assert abs(chi.root_number()) == pytest.approx(1)


def test_root_number_matches_lmfdb():
    assert conrey(5, 2).root_number() == pytest.approx(np.exp(2j * np.pi * 0.0881040955873917))


def test_root_number_needs_primitive():
    with pytest.raises(ValueError):
        conrey(12, 5).root_number()


@pytest.mark.parametrize("p", [3, 5, 7, 11, 13, 101, 997])
def test_kronecker_is_legendre_for_odd_primes(p):
    for a in range(-30, 30):
        euler = pow(a % p, (p - 1) // 2, p)
        assert kronecker(a, p) == (-1 if euler == p - 1 else euler)


def test_kronecker_is_multiplicative_in_n():
    for a in (-23, -4, 5, 12, 13):
        for m in range(1, 40):
            for n in range(1, 40):
                assert kronecker(a, m * n) == kronecker(a, m) * kronecker(a, n)


def test_fundamental_discriminants():
    neg = [d for d in range(-1, -30, -1) if is_fundamental(d)]
    pos = [d for d in range(30) if is_fundamental(d)]
    assert neg == [-3, -4, -7, -8, -11, -15, -19, -20, -23, -24]
    assert pos == [5, 8, 12, 13, 17, 21, 24, 28, 29]


@pytest.mark.parametrize("d", [-3, -4, -8, -23, 5, 8, 12, 229, -163])
def test_quadratic_character(d):
    chi = quadratic(d)
    assert chi.is_primitive and chi.is_real and chi.order == 2
    assert chi.parity == (d < 0)
    assert chi.root_number() == pytest.approx(1)


@pytest.mark.parametrize(
    "d, h", [(-3, 1), (-4, 1), (-23, 3), (-47, 5), (-56, 4), (-71, 7), (-163, 1), (-260, 8)]
)
def test_class_numbers(d, h):
    assert class_number(d) == h


def test_class_number_one_discriminants():
    heegner = [d for d in range(-3, -200, -1) if is_fundamental(d) and class_number(d) == 1]
    assert heegner == [-3, -4, -7, -8, -11, -19, -43, -67, -163]


@pytest.mark.parametrize(
    "d, unit", [(5, (1, 1)), (8, (2, 1)), (13, (3, 1)), (229, (15, 1)), (316, (160, 9)),
                (376, (4286590, 221064))],
)
def test_fundamental_units(d, unit):
    x, y = fundamental_unit(d)
    assert (x, y) == unit
    assert abs(x * x - d * y * y) == 4
