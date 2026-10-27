import json
from pathlib import Path

import mpmath

from zeta_explorer.explicit import psi_from_zeros
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
}
rows = ",\n".join(f" {json.dumps(k)}: {json.dumps(v)}" for k, v in out.items())
Path(__file__).with_name("reference.json").write_text("{\n" + rows + "\n}\n")
