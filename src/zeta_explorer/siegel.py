from __future__ import annotations

from dataclasses import dataclass

import mpmath
import numpy as np

CONVEXITY = {2: "2.97655", 4: "5.13781", 8: "9.49562"}
J0 = "1.912"


@dataclass(frozen=True)
class Constants:
    r: float
    phi: float
    E: float


def constants(r: float, dps: int = 30) -> Constants:
    with mpmath.workdps(dps):
        mpf, pi, log = mpmath.mpf, mpmath.pi, mpmath.log
        r = mpf(r)
        R = r + mpf(7) / 8
        C = {k: mpf(v) for k, v in CONVEXITY.items()}
        th = [pi / 2] + [mpmath.acos((r + mpf(o) / 8) / R) for o in (1, 2, 4, 6)] + [mpf(0)]
        s = [mpmath.sin(x) for x in th]
        c = [mpmath.cos(x) for x in th]

        def A(n, m):
            return s[m] - s[n]

        def B(n, m):
            return (th[n] - th[m]) / 2 + (s[n] * c[n] - s[m] * c[m]) / 2

        def weight(shift, n, m):
            return (r + shift) / (pi * R) * A(n, m) + B(n, m) / pi

        phi = r * A(1, 5) / (pi * R) + B(0, 1) / (pi * (8 * r + 1)) + B(1, 5) / pi
        rho = -2 * (A(0, 1) / (pi * R) + B(0, 1) / (pi * (r + mpf(1) / 8)))
        K = (2 * log(C[8]) + log(2 * r + 2) / 8) * B(0, 1) / (pi * (r + mpf(1) / 8))
        arc = abs(2 + r - R * mpmath.expjpi(th[3] / pi))
        K += 2 * (log(arc) - log(2 * pi)) * weight(mpf(1) / 2, 3, 5)
        for k in (2, 3):
            a = mpf(2) ** (k - 4)
            u = 2 ** (6 - k) * log(C[2 ** (4 - k)]) + 2 * log(mpf(15) / 8 + a + 2 * r)
            v = 2 ** (6 - k) * log(C[2 ** (5 - k)]) + log(mpf(15) / 8 + a + 2 * r)
            K += u * weight(a / 2, k - 1, k) - v * weight(a, k - 1, k)
        for k in (4, 5):
            a = mpf(2) ** (3 - k)
            u = 2 ** (k - 1) * log(C[2 ** (k - 3)]) + 2 * log(mpf(7) / 8 + a)
            v = 2 ** (k - 1) * log(C[2 ** (k - 2)]) + log(mpf(7) / 8 + a)
            K += v * weight(1 - a, k - 1, k) - u * weight(1 - a / 2, k - 1, k)
        K += mpf(J0) / (pi * R) - log(2 * pi) + mpmath.euler / 2 + 1
        K += mpmath.digamma((3 + r) / 2) / 2
        E = K + rho * log(1 + 1 / r)
        return Constants(float(r), float(phi), float(E))


def bound(q: np.ndarray | int, k: Constants, c: float = 0.2) -> np.ndarray:
    L = np.log(np.asarray(q, dtype=float))
    r = k.r
    R = r + 7 / 8
    return c / (r * (r * L + c)) + (r + c / L) / R**2 + k.phi * L + k.E
