"""Stack the RA maps of several processed sessions into one deeper map.

Each session's ra_map.npz (made by hi-process) holds the RA-binned spectra
and the number of cycles behind every bin. Stacking is a cycle-weighted
average per (RA bin, velocity channel), so N nights at the same elevation
give a map ~sqrt(N) deeper, and sessions covering different RA ranges (or
interrupted ones) simply fill in each other's gaps.

Only stack sessions taken at the SAME dish pointing -- mixing declinations
would average unrelated sky. A warning is printed if the per-bin
declinations disagree by more than a beam's fraction.

Usage:
    hi-stack data/raw/night1 data/raw/night2 ... [--out data/stacks/week1]
"""

from __future__ import annotations

import argparse
import pathlib
import time

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import extract, velocity
from .config import load_config
from .coverage import coverage_plot, find_strip_cache
from .process import _extraction_plots


def _session_dec(sdir: pathlib.Path):
    p = sdir / "products" / "ra_map.npz"
    if not p.exists():
        return None
    z = np.load(p)
    dec, counts = z["dec"], z["counts"]
    if (counts > 0).sum() == 0:
        return None
    return float(np.nanmean(dec[counts > 0]))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Stack processed sessions into a deeper RA map")
    ap.add_argument("sessions", nargs="*", help="session dirs (each must be hi-processed)")
    ap.add_argument("--auto", action="store_true",
                    help="discover all processed sessions under paths.raw_dir, "
                         "group them by declination (same dish setting), and "
                         "stack each group into data/stacks/dec<+XX>/")
    ap.add_argument("--dec-tol", type=float, default=2.0,
                    help="max dec difference [deg] for two sessions to count "
                         "as the same pointing in --auto mode")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--out", default=None,
                    help="output directory (default data/stacks/<timestamp>)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    if args.auto:
        raw = pathlib.Path(cfg.get("paths", {}).get("raw_dir", "data/raw"))
        cands = [(d, _session_dec(d)) for d in sorted(raw.iterdir())
                 if d.is_dir() and "sunscan" not in d.name.lower()] \
            if raw.exists() else []
        cands = [(d, dec) for d, dec in cands if dec is not None]
        if not cands:
            raise SystemExit("[stack] --auto found no processed sessions "
                             "(run hi-process --all first)")
        groups: list[list] = []
        for d, dec in cands:
            for g in groups:
                if abs(g[0][1] - dec) <= args.dec_tol:
                    g.append((d, dec))
                    break
            else:
                groups.append([(d, dec)])
        rc = 0
        for g in groups:
            dec0 = float(np.mean([dec for _, dec in g]))
            out = pathlib.Path("data/stacks") / f"dec{dec0:+05.1f}"
            print(f"[stack] group dec {dec0:+.1f}: "
                  f"{', '.join(d.name for d, _ in g)} -> {out}")
            rc |= _stack([str(d) for d, _ in g], cfg, out)
        return rc

    if not args.sessions:
        ap.error("no sessions given (pass directories or use --auto)")
    out = pathlib.Path(args.out) if args.out else \
        pathlib.Path("data/stacks") / time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    return _stack(args.sessions, cfg, out)


def _stack(sessions: list[str], cfg: dict, out: pathlib.Path) -> int:
    out.mkdir(parents=True, exist_ok=True)

    ra = v = None
    num = den = None
    dec_num = dec_den = None
    for sdir in sessions:
        p = pathlib.Path(sdir) / "products" / "ra_map.npz"
        if not p.exists():
            raise SystemExit(f"[stack] {p} missing -- run hi-process on {sdir} first")
        z = np.load(p)
        if ra is None:
            ra, v = z["ra"], z["v"]
            num = np.zeros((len(ra), len(v)))
            den = np.zeros((len(ra), len(v)))
            dec_num = np.zeros(len(ra))
            dec_den = np.zeros(len(ra))
        if len(z["ra"]) != len(ra) or not np.allclose(z["v"], v):
            raise SystemExit(f"[stack] {sdir} has a different RA/velocity grid -- "
                             "reprocess all sessions with the same config")
        W, counts, dec = z["W"], z["counts"], z["dec"]
        fin = np.isfinite(W)
        num += np.where(fin, W, 0.0) * counts[:, None]
        den += fin * counts[:, None]
        good = counts > 0
        dec_num[good] += dec[good] * counts[good]
        dec_den[good] += counts[good]
        print(f"[stack] {sdir}: {int(counts.sum())} cycles over "
              f"{int(good.sum())} RA bins")

    with np.errstate(invalid="ignore", divide="ignore"):
        W = (num / den).astype(np.float32)
        dec_b = dec_num / dec_den
    W[den == 0] = np.nan
    filled = dec_den > 0
    if filled.any():
        spread = np.nanmax(dec_b[filled]) - np.nanmin(dec_b[filled])
        if spread > 3.0:
            print(f"[stack] WARNING: declinations span {spread:.1f} deg -- "
                  "are these really the same pointing?")

    l_b = np.full(len(ra), np.nan)
    b_b = np.full(len(ra), np.nan)
    l_b[filled], b_b[filled] = velocity.galactic_lb(ra[filled], dec_b[filled])
    counts_total = dec_den.astype(int)
    np.savez_compressed(out / "ra_map.npz", ra=ra, dec=dec_b, l=l_b, b=b_b,
                        v=v, W=W, counts=counts_total)

    fig, ax = plt.subplots(figsize=(11, 5))
    vm = np.nanpercentile(W, [2, 99.5])
    im = ax.pcolormesh(ra, v, W.T, vmin=vm[0], vmax=vm[1], cmap="inferno",
                       shading="nearest")
    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("$v_{LSR}$ [km/s]")
    ax.set_title(f"stacked RA x velocity map ({len(sessions)} sessions, "
                 f"{int(counts_total.sum())} cycles)")
    fig.colorbar(im, ax=ax, label="$T_B$ [K]")
    fig.tight_layout()
    fig.savefig(out / "waterfall_ra.png", dpi=130)
    plt.close(fig)

    coverage_plot(out / "coverage.png", ra, dec_b, filled,
                  float(cfg["pointing"]["beam_fwhm_deg"]),
                  strip_path=find_strip_cache(cfg),
                  title=f"stacked coverage ({len(sessions)} sessions)")
    rows = extract.analyze_bins(v, W[filled], ra[filled], dec_b[filled],
                                l_b[filled], b_b[filled], out / "extraction.csv")
    _extraction_plots(out, v, W[filled], ra[filled], l_b[filled], b_b[filled], rows)
    print(f"[stack] stacked products in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
