from zeta_explorer.gf import factor, monics, squarefree
from zeta_explorer.hyperelliptic import irreducible_counts
from zeta_explorer.hypermoments import (
    irreducibles,
    squarefrees,
    squares,
)


def test_irreducibles():
    for q in (3, 5):
        for d in range(1, 5):
            Ps = irreducibles(q, d)
            assert len(Ps) == irreducible_counts(q, d)[d]
            assert all(P[-1] == 1 and len(P) == d + 1 for P in Ps)
            assert all(len(factor(P, q)) == 1 for P in Ps[:30])


def test_squarefrees_with_factor_degrees():
    q, n = 3, 5
    got = {f: degs for f, degs in squarefrees(q, n)}
    want = [f for d in range(n + 1) for f in monics(q, d) if squarefree(f, q)]
    assert set(got) == set(want) and len(got) == len(want)
    for f in want[1:200]:
        assert sorted(got[f]) == sorted(len(P) - 1 for P in factor(f, q))


def test_every_monic_is_squarefree_times_square():
    q = 3
    counts = [0] * 7
    hs = squares(q, 6)
    for f, _ in squarefrees(q, 6):
        for h, _ in hs:
            j = len(f) - 1 + 2 * (len(h) - 1)
            if j <= 6:
                counts[j] += 1
    assert counts == [q**j for j in range(7)]
