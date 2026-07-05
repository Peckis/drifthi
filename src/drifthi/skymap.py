"""All-sky coverage map: everything you have ever observed, on one plot.

Scans all processed sessions (and any explicitly given directories), paints
each drift track's beam smear onto a fine sky grid using exact spherical
caps, and draws the result in an Aitoff projection with the Milky Way band
(|b| < 10 deg) outlined. Reports the solid-angle fraction of the whole sky
and of the Milky Way band that your beam has covered.

Usage:
    hi-skymap                        # all processed sessions in paths.raw_dir
    hi-skymap data/raw/night1 ...    # only specific sessions
    hi-skymap --band-b 15 --out mymap.png
"""

from __future__ import annotations

import argparse
import pathlib

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .config import load_config

GRID_STEP = 0.5  # deg


def _tracks(session_dirs) -> list[tuple[str, np.ndarray, np.ndarray]]:
    out = []
    for d in session_dirs:
        p = pathlib.Path(d) / "products" / "ra_map.npz"
        if not p.exists():
            print(f"[skymap] {d}: not processed yet, skipping")
            continue
        z = np.load(p)
        sel = z["counts"] > 0
        if sel.any():
            out.append((pathlib.Path(d).name, z["ra"][sel], z["dec"][sel]))
    return out


def paint_coverage(tracks, beam_fwhm: float):
    """Boolean sky grid (nra, ndec) marking pixels within beam/2 of any track."""
    ra_g = np.arange(0.0, 360.0, GRID_STEP) + GRID_STEP / 2
    dec_g = np.arange(-90.0, 90.0, GRID_STEP) + GRID_STEP / 2
    cov = np.zeros((len(ra_g), len(dec_g)), dtype=bool)
    r = np.radians(beam_fwhm / 2.0)
    cosr = np.cos(r)
    for _, ras, decs in tracks:
        for ra0, dec0 in zip(ras, decs):
            j = (np.abs(dec_g - dec0) <= np.degrees(r))
            if not j.any():
                continue
            sd0, cd0 = np.sin(np.radians(dec0)), np.cos(np.radians(dec0))
            sd = np.sin(np.radians(dec_g[j]))
            cd = np.cos(np.radians(dec_g[j]))
            # spherical cap: half-width in RA at each dec row
            with np.errstate(invalid="ignore", divide="ignore"):
                cos_dra = np.clip((cosr - sd0 * sd) / np.maximum(cd0 * cd, 1e-9),
                                  -1.0, 1.0)
            dra = np.degrees(np.arccos(cos_dra))
            dist = np.abs((ra_g[:, None] - ra0 + 180.0) % 360.0 - 180.0)
            cov[:, j] |= dist <= dra[None, :]
    return ra_g, dec_g, cov


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="All-sky Aitoff coverage map")
    ap.add_argument("sessions", nargs="*",
                    help="session dirs (default: all processed under paths.raw_dir)")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--out", default="data/skymap.png")
    ap.add_argument("--band-b", type=float, default=10.0,
                    help="|b| defining the 'Milky Way band' for the statistics")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    beam = float(cfg["pointing"]["beam_fwhm_deg"])

    if args.sessions:
        dirs = args.sessions
    else:
        raw = pathlib.Path(cfg.get("paths", {}).get("raw_dir", "data/raw"))
        dirs = [d for d in sorted(raw.iterdir())
                if d.is_dir() and "sunscan" not in d.name.lower()] \
            if raw.exists() else []
    tracks = _tracks(dirs)
    if not tracks:
        raise SystemExit("[skymap] no processed sessions found -- run hi-process first")
    print(f"[skymap] {len(tracks)} sessions, beam FWHM {beam:.0f} deg")

    ra_g, dec_g, cov = paint_coverage(tracks, beam)

    # galactic latitude of every grid pixel (for the Milky Way band)
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    rr, dd = np.meshgrid(ra_g, dec_g, indexing="ij")
    b = SkyCoord(rr.ravel() * u.deg, dd.ravel() * u.deg).galactic.b.deg
    b = b.reshape(rr.shape)
    mw = np.abs(b) < args.band_b

    w = np.cos(np.radians(dd))                       # solid angle weight
    sky_pct = 100.0 * np.sum(cov * w) / np.sum(w)
    mw_pct = 100.0 * np.sum(cov * mw * w) / max(np.sum(mw * w), 1e-9)
    print(f"[skymap] sky covered: {sky_pct:.1f}%   "
          f"Milky Way band (|b|<{args.band_b:.0f}) covered: {mw_pct:.1f}%")

    # ---- Aitoff plot: x = 180 - RA (astronomical convention, RA grows left)
    x = np.radians(180.0 - ra_g)
    order = np.argsort(x)
    X, Y = np.meshgrid(x[order], np.radians(dec_g), indexing="ij")

    fig = plt.figure(figsize=(13, 7))
    ax = fig.add_subplot(111, projection="aitoff")
    ax.grid(alpha=0.3)
    ax.pcolormesh(X, Y, np.ma.masked_where(~mw[order], np.ones_like(X)),
                  cmap=matplotlib.colors.ListedColormap(["#d8c9a3"]),
                  alpha=0.5, shading="nearest", rasterized=True)
    ax.pcolormesh(X, Y, np.ma.masked_where(~cov[order], np.ones_like(X)),
                  cmap=matplotlib.colors.ListedColormap(["#00b7c7"]),
                  alpha=0.8, shading="nearest", rasterized=True)

    # galactic plane line and center
    gal = SkyCoord(l=np.linspace(0, 360, 720) * u.deg, b=0 * u.deg,
                   frame="galactic").icrs
    gx = np.radians(180.0 - gal.ra.deg)
    gy = np.radians(gal.dec.deg)
    cut = np.where(np.abs(np.diff(gx)) > np.pi / 2)[0] + 1  # split at wrap
    for seg_x, seg_y in zip(np.split(gx, cut), np.split(gy, cut)):
        ax.plot(seg_x, seg_y, "k--", lw=0.9)
    gc = SkyCoord(l=0 * u.deg, b=0 * u.deg, frame="galactic").icrs
    ax.plot(np.radians(180 - gc.ra.deg), np.radians(gc.dec.deg), "k*", ms=12)

    ax.set_xticklabels([f"{(180 - d) % 360:.0f}°"
                        for d in np.arange(-150, 180, 30)], fontsize=8)
    ax.set_title(f"observed sky (beam-smeared tracks, {len(tracks)} sessions)\n"
                 f"sky: {sky_pct:.1f}%   Milky Way band |b|<{args.band_b:.0f}°: "
                 f"{mw_pct:.1f}%   (RA/Dec, equatorial)", fontsize=11)

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"[skymap] wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
