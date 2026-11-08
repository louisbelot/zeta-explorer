import json
import math

import mpmath
import numpy as np
import pytest
from flint import arb

from zeta_explorer import certified as C
from zeta_explorer.riemann_siegel import theta, z_exact, z_rs
from zeta_explorer.zeros import gram_points
from zeta_explorer.zeta import zeta


@pytest.mark.parametrize("t", [0.0, 3.5, 14.1, 100.0, 1000.0, 12345.6])
def test_z_ball_contains_mpmath(t):
    z = C.z_ball(t, 128)
    with mpmath.workdps(40):
        ref = mpmath.siegelz(mpmath.mpf(t))
    assert z.rad() < 1e-30
    assert abs(float(z.mid()) - float(ref)) < 1e-14 * max(1.0, abs(float(ref)))


@pytest.mark.parametrize("t, tol", [(50.0, 1e-12), (250.0, 5e-9), (1e3, 1e-10), (1e4, 1e-10),
                                    (1e5, 1e-9)])
def test_own_z_against_ball(t, tol):
    ts = t + np.linspace(0, 3, 13)
    own = np.array([z_exact(x) for x in ts]) if t < 200 else z_rs(ts)
    for x, z in zip(ts, own, strict=True):
        assert abs(z - float(C.z_ball(float(x), 128).mid())) < tol


@pytest.mark.parametrize("s", [2 + 0j, 0.5 + 30j, -3.5 + 7j, 0.9 + 500j, 3 - 40j])
def test_own_zeta_inside_ball(s):
    ball = C.zeta_ball(s, 128)
    ref = complex(float(ball.real.mid()), float(ball.imag.mid()))
    assert abs(zeta(s) - ref) < 1e-12 * max(1.0, abs(ref))


def test_theta_ball_matches_float():
    ts = np.array([10.0, 100.0, 1e3, 1e5])
    for t, th in zip(ts, theta(ts), strict=True):
        assert abs(float(C.theta_ball(float(t), 128).mid()) - th) < 1e-12 * max(1.0, abs(th))


def test_z_sign_splits_float_neighbours_of_first_zero():
    assert C.z_sign(14.134725141734693) == -1
    assert C.z_sign(14.134725141734694) == 1


@pytest.mark.parametrize("n", [-1, 0, 1, 100, 1000, 20000])
def test_gram_ball_contains_arb_gram_point(n):
    guess = C.G_MINUS_1 if n < 0 else float(gram_points([n])[0])
    g = C.gram_ball(n, guess)
    assert g.rad() < 1e-30
    with C.precision(128):
        if n >= 0:
            assert g.overlaps(arb.gram_point(n))
        assert abs(C.theta_ball(g.mid(), 128) - n * arb.pi()) < 1e-30


def test_sign_changes_skip_undecided_points():
    assert C.count_sign_changes([1, 0, -1, -1, 0, 1]) == 2
    assert C.count_sign_changes([1, 1, 1]) == 0


def test_blocks_needed_from_trudgian():
    for g in (600.0, 1e4, 1e8):
        L = math.log(g)
        assert C.blocks_needed(g) >= L * (2.067 + 0.059 * L) / (6 * math.pi)
        assert C.blocks_needed(g) - 1 < L * (2.067 + 0.059 * L) / (6 * math.pi)


@pytest.mark.parametrize("T", [100.0, 1000.0, 3000.0])
def test_certify_matches_arb_count(T):
    c = C.certify(T)
    assert float(c.g_n.mid()) >= max(T, C.T_MIN)
    assert c.zeros == c.n + 1
    with C.precision(128):
        assert int(c.g_n.mid().zeta_nzeros().unique_fmpz()) == c.zeros
    assert c.turing_blocks >= C.blocks_needed(float(c.g_n.mid()))


def test_certificate_log_is_consistent():
    c = C.certify(600.0)
    log = json.loads(json.dumps(c.log()))
    assert log["zeros_on_line"] == c.n + 1
    blocks = log["blocks"]
    assert blocks[0]["j"] == -1
    assert all(a["k"] == b["j"] for a, b in zip(blocks[:-1], blocks[1:], strict=True))
    assert blocks[-1]["k"] == c.n
    assert sum(C.count_sign_changes(b["signs"]) for b in blocks) == c.n + 1
    for b in blocks:
        assert len(b["points"]) == len(b["signs"])
        assert b["points"] == sorted(b["points"])
