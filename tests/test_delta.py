import math

import mpmath as mp
import numpy as np
import pytest

from zeta_explorer import delta as D
from zeta_explorer.explicit import primes

LMFDB_TAU = [1, -24, 252, -1472, 4830, -6048, -16744, 84480, -113643, -115920, 534612,
             -370944, -577738, 401856, 1217160, 987136, -6905934]
LMFDB_ZEROS = [9.222379399921102522, 13.907549861392134406, 17.442776978234473314,
               19.656513141954961000, 22.336103637209867276, 25.274636548112365354,
               26.804391158350403033, 28.831682624186875445, 31.178209498360259064,
               32.774875382231207442]


def test_tau_lmfdb():
    assert list(D.tau(17)[1:]) == LMFDB_TAU


def test_tau_against_flint_product():
    flint = pytest.importorskip("flint")
    N = 1500
    eta = flint.fmpz_poly([1])
    for n in range(1, N):
        eta = (eta * flint.fmpz_poly([1] + [0] * (n - 1) + [-1])).truncate(N)
    ref = [int(c) for c in eta.pow_trunc(24, N).coeffs()]
    assert list(D.tau(N)[1:]) == ref


def test_ramanujan_691():
    N = 3000
    t = D.tau(N)
    sigma = [0] * (N + 1)
    for d in range(1, N + 1):
        for m in range(d, N + 1, d):
            sigma[m] += pow(d, 11, 691)
    assert all((t[n] - sigma[n]) % 691 == 0 for n in range(1, N + 1))


def test_hecke_relations():
    t = D.tau(2000)
    for m in range(2, 40):
        for n in range(2, 2000 // m + 1):
            if math.gcd(m, n) == 1:
                assert t[m * n] == t[m] * t[n]
    for p in primes(44):
        assert t[p * p] == t[p] ** 2 - p**11


def test_deligne_bound_exact():
    t = D.tau(1 << 15)
    for p in primes(1 << 15):
        assert t[p] ** 2 < 4 * int(p) ** 11


def test_moduli_cover_the_bound():
    for N in (1024, 1 << 20):
        assert math.prod(D._moduli(N)) > 4 * N**6


def test_lam_is_normalized_tau():
    n = np.arange(1, 101)
    assert np.allclose(D.lam(100)[1:] * n**5.5, [float(x) for x in D.tau(100)[1:]], rtol=1e-15)


def test_zeros_lmfdb():
    assert np.allclose(D.zeros(33.5), LMFDB_ZEROS, rtol=1e-13)


def test_zero_count_follows_theta():
    z = D.zeros(200)
    for T in (100, 200):
        n = int(np.sum(z < T))
        assert abs(n - D.theta(T) / math.pi) < 1
    assert int(np.sum(z < 200)) == 159


def test_central_value_lmfdb():
    assert D.l_delta(0.5).real == pytest.approx(0.7921228386, abs=1e-10)


def _z_mpmath(t, dps):
    with mp.workdps(dps):
        s = mp.mpf(0.5) + 1j * mp.mpf(t)
        A = 1 / (2 * mp.pi)
        tot = mp.mpc(0)
        for n, tn in enumerate(D.tau(140)[1:], start=1):
            tot += mp.mpf(int(tn)) / mp.mpf(n) ** 5.5 * (A / n) ** s * mp.gammainc(s + 5.5, n / A)
        return float(2 * mp.re(tot) / (mp.sqrt(A) * abs(mp.gamma(6 + 1j * mp.mpf(t)))))


@pytest.mark.parametrize("t, dps", [(5.0, 40), (50.0, 70), (100.0, 110)])
def test_hardy_z_against_mpmath(t, dps):
    assert D.hardy_z(t) == pytest.approx(_z_mpmath(t, dps), abs=1e-12)


def test_split_point_and_functional_equation():
    s = 0.3 + 7j
    assert D.completed(s) == pytest.approx(D.completed(s, 1.3 + 0.4j), rel=1e-12)
    assert D.completed(s) == pytest.approx(D.completed(1 - s), rel=1e-12)


def test_period_ratios():
    r = [math.comb(10, n) * D.completed(n - 4.5).real for n in range(11)]
    even = np.array(r[0::2]) / r[2]
    odd = np.array(r[1::2]) / r[1] * 4
    assert np.allclose(even, [36 / 691, 1, 3, 3, 1, 36 / 691], rtol=1e-12)
    assert np.allclose(odd, [4, 25, 42, 25, 4], rtol=1e-12)


def test_completed_strip_check():
    with pytest.raises(ValueError):
        D.completed(7.0)


def test_sato_tate_small():
    st = D.sato_tate(1 << 15)
    assert np.all(np.abs(st.x) < 1)
    assert st.moment(2) == pytest.approx(D.catalan(1), abs=0.05)
    assert st.moment(4) == pytest.approx(D.catalan(2), abs=0.15)
    assert abs(st.moment(1)) < 0.05
    assert st.ks() * math.sqrt(len(st.p)) < 1.36


def test_st_cdf_and_density():
    th = np.linspace(0.1, 3.0, 9)
    h = 1e-6
    assert np.allclose((D.st_cdf(th + h) - D.st_cdf(th - h)) / (2 * h), D.st_density(th))
    assert D.st_cdf(np.pi) == pytest.approx(1.0)
    assert [D.catalan(k) for k in range(6)] == [1, 1, 2, 5, 14, 42]


def test_psi_against_prime_powers():
    X = 5000
    st = D.sato_tate(X)
    lam = D.lam(X)
    total = 0.0
    for p in primes(X):
        s0, s1 = 2.0, lam[p]
        q = p
        while q <= X:
            total += s1 * math.log(p)
            s0, s1 = s1, lam[p] * s1 - s0
            q *= p
    assert st.psi(1, X + 0.5) == pytest.approx(total, rel=1e-12)


def test_root_numbers():
    assert [D.root_number(m) for m in range(1, 9)] == [1, 1, -1, 1, -1, 1, 1, 1]
