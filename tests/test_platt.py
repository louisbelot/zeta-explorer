from pathlib import Path

import pytest

pytest.importorskip("flint")

from zeta_explorer import platt  # noqa: E402

DATA = Path(__file__).parent / "data"


def block(e: int) -> str:
    return (DATA / f"platt_1e{e}.txt").read_text()


def test_parse():
    ns, gs = platt.parse(block(9))
    assert ns[0] == 10**9 and ns[-1] == 10**9 + 49
    assert gs[0] == "371870203.8370280527340547959866251910007"


@pytest.mark.parametrize("e", [9, 10, 11])
def test_platt_block(e):
    c = platt.check(block(e))
    assert c.all_isolated
    assert c.max_diff < 1e-30
    assert c.window.complete
    assert c.window.a < 10**e - 1 and c.window.b >= 10**e + 48
    assert c.indices_match


def test_wrong_digit_is_not_isolated():
    _, gs = platt.parse(block(9))
    bad = gs[3][:-12] + "9" + gs[3][-11:]
    assert platt.isolate(gs[3]) is not None
    assert platt.isolate(bad, delta="1e-28") is None


def test_missing_zero_breaks_indices():
    ns, gs = platt.parse(block(9))
    gs = gs[:20] + gs[21:]
    c = platt.check("\n".join(f"{n} {g}" for n, g in zip(ns, gs, strict=False)))
    assert c.window.complete
    assert not c.indices_match
