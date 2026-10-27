import json
from pathlib import Path

import mpmath
import numpy as np
import pytest

from zeta_explorer.riemann_siegel import rs_coefficients
from zeta_explorer.site import build, rs_data, spacing_data, zeros_data

WEB_DATA = Path(__file__).parents[1] / "web" / "data"


def _horner(c, z):
    y = np.zeros_like(z)
    for a in reversed(c):
        y = y * z + a
    return y


def test_rs_data_truncation_is_harmless():
    z = np.linspace(-0.5, 0.5, 2001)
    for short, full in zip(rs_data(), rs_coefficients(), strict=True):
        assert len(short) < len(full)
        assert np.max(np.abs(_horner(short, z) - np.polyval(full, z))) < 1e-16


def test_zeros_data_matches_mpmath():
    gamma = zeros_data(40)["gamma"]
    assert len(gamma) == 40
    for n in (1, 2, 17, 40):
        assert gamma[n - 1] == pytest.approx(float(mpmath.zetazero(n).imag), abs=1e-9)


def test_spacing_data_is_normalized():
    d = spacing_data(10**6, 2000, bins=35)
    assert d["count"] == 2000
    width = np.diff(d["spacing"]["edges"])
    mass = np.sum(np.array(d["spacing"]["density"]) * width)
    assert 0.99 < mass <= 1.0 + 1e-6
    assert len(d["pair"]["u"]) == len(d["pair"]["r"]) == 35
    assert 0 < d["ks"] < 0.05


def test_build_writes_json(tmp_path):
    names = build(tmp_path, n_zeros=30, spacing_n0=10**6, spacing_count=500)
    assert names == ["rs.json", "spacing.json", "zeros.json"]
    assert json.loads((tmp_path / "zeros.json").read_text())["count"] == 30


def test_committed_data_is_current():
    assert json.loads((WEB_DATA / "rs.json").read_text()) == rs_data()
    zeros = json.loads((WEB_DATA / "zeros.json").read_text())
    assert zeros["gamma"][:200] == zeros_data(200)["gamma"]
