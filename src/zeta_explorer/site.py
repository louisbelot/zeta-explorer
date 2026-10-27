from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from zeta_explorer.multieval import zeros_fast
from zeta_explorer.riemann_siegel import rs_coefficients
from zeta_explorer.spacing import (
    first_zeros,
    gue_spacing,
    gue_spacing_cdf,
    ks_distance,
    montgomery,
    pair_correlation,
    wigner_surmise,
)

RS_TOL = 1e-18


def rs_data(tol: float = RS_TOL) -> list[list[float]]:
    out = []
    for c in rs_coefficients():
        ascending = c[::-1]
        size = np.abs(ascending) * 0.5 ** np.arange(len(ascending))
        keep = np.flatnonzero(size > tol)
        out.append(ascending[: keep[-1] + 1].tolist() if len(keep) else [0.0])
    return out


def zeros_data(n: int) -> dict:
    gamma = first_zeros(n)
    return {"count": n, "gamma": [round(float(g), 10) for g in gamma]}


def spacing_data(n0: int, count: int, bins: int = 50, s_max: float = 3.5) -> dict:
    win, u = zeros_fast(n0, n0 + count)
    x = win.unfold(u)
    gaps = np.diff(x)
    hist, edges = np.histogram(gaps, bins=bins, range=(0, s_max), density=True)
    s = np.linspace(0, s_max, 141)
    r_mid, r = pair_correlation(x, u_max=3.0, bins=bins)
    uu = np.linspace(0, 3.0, 121)
    return {
        "n0": n0,
        "count": len(u),
        "t0": win.t0,
        "ks": ks_distance(gaps, gue_spacing_cdf),
        "spacing": {
            "edges": edges.round(6).tolist(),
            "density": hist.round(6).tolist(),
            "s": s.round(6).tolist(),
            "gue": gue_spacing(s).round(6).tolist(),
            "wigner": wigner_surmise(s).round(6).tolist(),
            "poisson": np.exp(-s).round(6).tolist(),
        },
        "pair": {
            "u": r_mid.round(6).tolist(),
            "r": r.round(6).tolist(),
            "uu": uu.round(6).tolist(),
            "montgomery": montgomery(uu).round(6).tolist(),
        },
    }


def build(
    out: Path, n_zeros: int = 10_000, spacing_n0: int = 10**12, spacing_count: int = 40_000
) -> list[str]:
    out.mkdir(parents=True, exist_ok=True)
    files = {
        "rs.json": rs_data(),
        "zeros.json": zeros_data(n_zeros),
        "spacing.json": spacing_data(spacing_n0, spacing_count),
    }
    for name, data in files.items():
        (out / name).write_text(json.dumps(data, separators=(",", ":")))
    return sorted(files)


def main() -> None:
    p = argparse.ArgumentParser(description="write the data files for the web page")
    p.add_argument("out", type=Path, nargs="?", default=Path("web/data"))
    p.add_argument("--zeros", type=int, default=10_000)
    p.add_argument("--spacing-n0", type=int, default=10**12)
    p.add_argument("--spacing-count", type=int, default=40_000)
    args = p.parse_args()
    for name in build(args.out, args.zeros, args.spacing_n0, args.spacing_count):
        print(args.out / name)


if __name__ == "__main__":
    main()
