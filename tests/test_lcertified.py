import math

import numpy as np
import pytest
from flint import dirichlet_char

from zeta_explorer.characters import conrey, primitive_characters
from zeta_explorer.lcertified import (
    Certificate,
    bmor,
    certify,
    characters,
    main_term,
    refine,
    window_area,
    z_sign,
)
from zeta_explorer.lfunction import count_zeros_box, l_zeros


def test_characters_up_to_conjugation():
    for q in range(3, 40):
        prim = primitive_characters(q)
        real = sum(1 for _, c in prim if c.is_real)
        pairs = real + (len(prim) - real) // 2
        assert len(characters(q)) == pairs


def test_conrey_labels_agree():
    for q in (5, 8, 12, 15, 21, 27):
        for chi in characters(q):
            ours = conrey(q, int(chi.number()))
            for m in range(1, q):
                if math.gcd(m, q) == 1:
                    assert complex(chi(m)) == pytest.approx(ours(m), abs=1e-12)


@pytest.mark.parametrize("q,n,gamma", [(4, 3, 6.020948904697597), (3, 2, 8.039737155681)])
def test_first_zero(q, n, gamma):
    chi = dirichlet_char(q, n)
    a, b = refine(chi, gamma - 0.1, gamma + 0.1)
    assert (a + b) / 2 == pytest.approx(gamma, abs=1e-11)


def test_first_zero_against_float_finder():
    chi = dirichlet_char(5, 2)
    cert = certify(chi, 20.0)
    lo, hi = min((iv for iv in cert.zeros if iv[0] > 0), key=lambda iv: iv[0])
    a, b = refine(chi, lo, hi)
    floats = l_zeros(conrey(5, 2), 0.0, 10.0)
    assert (a + b) / 2 == pytest.approx(floats[floats > 0].min(), abs=1e-9)


def test_negative_side_is_the_conjugate():
    chi, bar = dirichlet_char(13, 2), dirichlet_char(13, 7)
    assert (2 * 7) % 13 == 1
    for t in (3.0, 7.5, 12.25):
        assert float(abs(chi.hardy_z(-t)).mid()) == pytest.approx(
            float(abs(bar.hardy_z(t)).mid()), rel=1e-12
        )


@pytest.mark.parametrize("a", [0, 1])
def test_main_term_near_smooth(a):
    for T in (5 / 7, 1.0, 3.0, 10.0, 50.0):
        smooth = bmor(T, 7, a)[0]
        assert abs(main_term(T, 7, a) - smooth) <= (2 - a) / (50 * T)


@pytest.mark.parametrize("q,n", [(4, 3), (5, 2), (8, 5), (12, 11), (13, 2)])
def test_count_against_argument_principle(q, n):
    cert = certify(dirichlet_char(q, n), 20.0)
    assert cert.complete
    T = next(T for T in (20.0, 19.9, 19.8) if not cert.straddles(T))
    assert cert.count(T) == count_zeros_box(conrey(q, n), T)


def test_counts_respect_bmor():
    for q in (3, 4, 7, 11, 24):
        for chi in characters(q):
            cert = certify(chi, 30.0)
            assert cert.complete
            for T in np.arange(1.0, 30.0, 0.5):
                if cert.straddles(float(T)):
                    continue
                N = cert.count(float(T))
                smooth, err, _ = bmor(float(T), q, cert.a)
                if math.log(q * (T + 2) / (2 * math.pi)) <= 1.567:
                    assert N == 0
                else:
                    assert abs(N - smooth) <= err


def test_dropped_zeros_show_up_as_missing():
    cert = certify(dirichlet_char(7, 3), 25.0)
    assert cert.complete
    kept = sorted(cert.zeros, key=lambda iv: abs(iv[0]))[2:]
    area = window_area(kept, cert.q, cert.a, cert.t1, cert.t2)
    worse = Certificate(cert.q, cert.n, cert.a, cert.t1, cert.t2, kept, area, cert.bound)
    assert worse.missing == pytest.approx(cert.missing + 2, abs=1e-9)
    assert not worse.complete


def test_z_sign_certified():
    chi = dirichlet_char(4, 3)
    assert z_sign(chi, 6.0) == 1 and z_sign(chi, 6.05) == -1
