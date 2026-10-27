import math

import numpy as np
import pytest

from zeta_explorer.characters import (
    characters,
    conrey,
    conrey_generator,
    factor,
    primitive_characters,
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
