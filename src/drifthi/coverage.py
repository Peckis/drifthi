"""Sky-coverage footprint plot.

Draws the patch of sky a session or stack actually covered: the beam swept
along the drift track produces an elongated band (not circles -- the sky
moves during the observation). If an HI4PI strip cache is available it is
drawn underneath as a column-density map, so you see your footprint on the
real sky.
"""

from __future__ import annotations

import pathlib

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def find_strip_cache(cfg: dict) -> pathlib.Path | None:
    d = pathlib.Path(cfg["hi4pi"]["cache_dir"])
    if not d.exists():
        return None
    caches = sorted(d.glob("*.npz"), key=lambda p: p.stat().st_size, reverse=True)
    return caches[0] if caches else None


def coverage_plot(out_png: pathlib.Path, ra_b, dec_b, filled, beam_fwhm: float,
                  strip_path: pathlib.Path | None = None, title: str = "sky coverage"):
    ra_f = np.asarray(ra_b)[filled]
    dec_f = np.asarray(dec_b)[filled]
    if ra_f.size == 0:
        return
    order = np.argsort(ra_f)
    ra_f, dec_f = ra_f[order], dec_f[order]
    half = beam_fwhm / 2.0

    fig, ax = plt.subplots(figsize=(11, 6))
    if strip_path is not None and pathlib.Path(strip_path).exists():
        from .hi4pi import load_strip
        ra_g, dec_g, v_g, cube = load_strip(strip_path)
        dv = float(np.median(np.diff(v_g)))
        nhi = 1.823e18 * np.nansum(np.nan_to_num(cube), axis=2) * dv
        im = ax.pcolormesh(ra_g, dec_g, np.log10(np.maximum(nhi.T, 1e19)),
                           cmap="magma", shading="nearest")
        fig.colorbar(im, ax=ax, label=r"HI4PI log10 $N_{HI}$ [cm$^{-2}$]")

    # the swept-beam band: beam smeared along the drift track
    ax.fill_between(ra_f, dec_f - half, dec_f + half, color="cyan", alpha=0.25)
    ax.plot(ra_f, dec_f - half, "cyan", lw=1.2)
    ax.plot(ra_f, dec_f + half, "cyan", lw=1.2)
    ax.plot(ra_f, dec_f, "cyan", lw=0.6, ls="--")
    # rounded end caps (RA radius stretched by 1/cos(dec))
    th = np.linspace(0, 2 * np.pi, 100)
    for re_, de_ in ((ra_f[0], dec_f[0]), (ra_f[-1], dec_f[-1])):
        ax.plot(re_ + half * np.cos(th) / max(np.cos(np.radians(de_)), 0.2),
                de_ + half * np.sin(th), "cyan", lw=1.2)

    ax.set_xlabel("RA [deg]")
    ax.set_ylabel("Dec [deg]")
    ax.invert_xaxis()
    ax.set_title(f"{title} (beam FWHM {beam_fwhm:.0f} deg swept along the track)")
    if strip_path is None or not pathlib.Path(strip_path).exists():
        ax.grid(alpha=0.3)
        ax.set_ylim(np.nanmin(dec_f) - beam_fwhm, np.nanmax(dec_f) + beam_fwhm)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
