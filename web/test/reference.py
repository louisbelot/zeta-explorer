import json
from pathlib import Path

import mpmath

from zeta_explorer.characters import conrey
from zeta_explorer.explicit import psi_from_zeros
from zeta_explorer.lfunction import l_z, verify_zeros
from zeta_explorer.spacing import first_zeros

mpmath.mp.dps = 40

ZETA_AT = [
    (2, 0), (0.5, 14.134725141734693), (0.5, 0), (0, 0), (-1, 0), (-2, 0), (-3.5, 7),
    (0.3, 50.2), (3, -20), (-10, 0.5), (0.999, 0.001), (0.5, 100), (1.5, 40), (-4.5, -30),
    (0.75, -3), (-0.5, 60), (4, 2), (0.01, 0.01),
]
Z_AT = [0, 1, 7.5, 14.134725141734693, 50, 123.4, 199.9, 200.1, 1000, 5000, 1e4, 1e5, 1e6]
THETA_AT = [0.5, 3, 10, 50, 200, 1e4, 1e7]
GRAM_AT = [0, 1, 10, 100, 1000, 100_000]
ZERO_RANGES = [(1, 29), (1990, 2000)]
PSI_AT = [2.5, 10.5, 30.5, 99.5]
CONREY_AT = [(3, 2), (4, 3), (5, 2), (7, 3), (8, 5), (9, 4), (12, 11), (16, 3), (24, 5), (25, 7),
             (40, 3), (63, 5), (64, 3), (99, 2)]
L_AT = [(3, 2, 0.5, 8.0397), (5, 2, 0.5, 30), (7, 3, 0.5, -12.5), (8, 5, 2, 1), (9, 4, 0.25, 70),
        (12, 11, 0.5, 0.1), (16, 3, 0.9, -45), (25, 2, 0.5, 3), (40, 3, 0.1, 22),
        (63, 5, 0.5, -60), (99, 2, 0.5, 15), (101, 2, 0.5, 150)]
LZ_AT = [(5, 2), (7, 3), (13, 2), (24, 5), (61, 2)]
VERIFY_AT = [(5, 2, 40), (7, 3, 30), (13, 2, 25), (40, 3, 20), (99, 2, 12)]


def _values(q, n):
    chi = conrey(q, n)
    return [[float(v.real), float(v.imag)] for v in chi.values]


def _dirichlet(q, n, s, t):
    z = mpmath.dirichlet(mpmath.mpc(s, t), [mpmath.mpc(*v) for v in _values(q, n)])
    return [q, n, s, t, float(z.real), float(z.imag)]


def _l_z(q, n):
    t = [-37.5, -8.2, 0.0, 4.4, 19.9, 64.0]
    return [q, n, t, [float(z) for z in l_z(t, conrey(q, n))]]


def _verify(q, n, T):
    v = verify_zeros(conrey(q, n), T)
    assert v.ok
    return [q, n, T, v.in_box, [float(z) for z in v.zeros]]


def _zeta(s, t):
    z = mpmath.zeta(mpmath.mpc(s, t))
    return [s, t, float(z.real), float(z.imag)]


out = {
    "zeta": [_zeta(s, t) for s, t in ZETA_AT],
    "Z": [[t, float(mpmath.siegelz(t))] for t in Z_AT],
    "theta": [[t, float(mpmath.siegeltheta(t))] for t in THETA_AT],
    "gram": [[n, float(mpmath.grampoint(n))] for n in GRAM_AT],
    "zeros": {f"{a}-{b}": [float(mpmath.zetazero(n).imag) for n in range(a, b + 1)]
              for a, b in ZERO_RANGES},
    "psi100": [[x, float(y)] for x, y in zip(PSI_AT, psi_from_zeros(PSI_AT, first_zeros(100)),
                                             strict=True)],
    "conrey": [[q, n, _values(q, n)] for q, n in CONREY_AT],
    "dirichlet": [_dirichlet(*a) for a in L_AT],
    "lz": [_l_z(q, n) for q, n in LZ_AT],
    "lzeros": [_verify(*a) for a in VERIFY_AT],
}
rows = ",\n".join(f" {json.dumps(k)}: {json.dumps(v)}" for k, v in out.items())
Path(__file__).with_name("reference.json").write_text("{\n" + rows + "\n}\n")
