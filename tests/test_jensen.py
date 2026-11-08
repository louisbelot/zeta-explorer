import mpmath

from zeta_explorer.jensen import taylor_gamma, theta_moments


def test_theta_moments_match_quadrature():
    ms = [0, 3, 40]
    F = theta_moments(ms, dps=40)
    with mpmath.workdps(50):
        theta = lambda t: mpmath.jtheta(3, 0, mpmath.exp(-mpmath.pi * t)) / 2 - 0.5  # noqa: E731
        pts = [mpmath.mpf(k) / 8 for k in range(81)]
        for m, f in zip(ms, F, strict=True):
            ref = mpmath.quad(lambda u, m=m: u**m * mpmath.exp(u / 4) * theta(mpmath.exp(u)), pts)
            assert abs(f / ref - 1) < 1e-35


def test_gamma_against_gorz_table():
    g = taylor_gamma([10, 100, 1000], dps=30)
    table = ["1.6323380490e-17", "6.5777263785e-205", "3.8760340890e-2567"]
    for x, ref in zip(g, table, strict=True):
        assert abs(x / mpmath.mpf(ref) - 1) < 1e-10


def test_gamma_against_taylor_series():
    g = taylor_gamma(range(4), dps=30)
    with mpmath.workdps(30):
        half = mpmath.mpf(1) / 2
        f = lambda z: (4 * z * z - 1) * mpmath.pi ** (-(half + z) / 2) * mpmath.gamma(  # noqa: E731
            (half + z) / 2) * mpmath.zeta(half + z)
        c = mpmath.taylor(f, 0, 6)
        for n in range(4):
            assert abs(g[n] / (c[2 * n] * mpmath.factorial(n)) - 1) < 1e-20
