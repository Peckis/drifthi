"""Compare an observed session with HI4PI and fit the true pointing.

Physics of the fit: a fixed az/el drift scan always traces a line of constant
declination, offset in hour angle from the meridian. Errors in (az, el) map
into (dec, HA) of the track, so cross-correlating the observed RA x velocity
waterfall against HI4PI (smoothed to the dish beam) over a grid of (daz, del)
offsets recovers the true pointing -- typically to a fraction of the beam,
because galactic-plane crossings and velocity structure are very distinctive.

The same fit yields an amplitude scale factor, i.e. an absolute intensity
calibration: T_sys,true = T_sys,assumed / gain.

Usage:  hi-compare data/raw/<session> [--config config.yaml]
Needs the strip cache from hi-fetch-hi4pi (or a toy cache from hi-simulate).
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import hi4pi, velocity
from .config import load_config


def _load_calibrated(session_dir: pathlib.Path):
    p = session_dir / "products" / "calibrated.npz"
    if not p.exists():
        raise FileNotFoundError(f"{p} not found -- run hi-process first")
    z = np.load(p, allow_pickle=False)
    return z["t"], z["v"], z["T"]


def _time_bin(t, T, bin_s: float):
    edges = np.arange(t[0], t[-1] + bin_s, bin_s)
    idx = np.digitize(t, edges) - 1
    tb, Ob = [], []
    for b in range(len(edges) - 1):
        sel = idx == b
        if sel.sum() == 0:
            continue
        tb.append(t[sel].mean())
        Ob.append(np.nanmean(T[sel], axis=0))
    return np.array(tb), np.array(Ob)


def _score(O, M):
    """(pearson r, gain) over pixels finite in both."""
    m = np.isfinite(O) & np.isfinite(M)
    if m.sum() < 200:
        return -np.inf, np.nan
    o, mm = O[m], M[m]
    o = o - o.mean()
    mc = mm - mm.mean()
    denom = np.sqrt((o**2).sum() * (mc**2).sum())
    if denom <= 0:
        return -np.inf, np.nan
    r = float((o * mc).sum() / denom)
    gain = float((O[m] * mm).sum() / max((mm**2).sum(), 1e-12))
    return r, gain


def fit_pointing(cfg, strip_path: pathlib.Path, t, v, T):
    pc = cfg["processing"]
    cc = cfg["compare"]
    tb, O = _time_bin(t, T, float(pc["time_bin_s"]))
    print(f"[compare] observed waterfall binned to {len(tb)} x {len(v)}")

    ra_g, dec_g, v_g, cube = hi4pi.load_strip(strip_path)
    sm = hi4pi.beam_smooth_strip(ra_g, dec_g, cube,
                                 float(cfg["pointing"]["beam_fwhm_deg"]))
    ev = hi4pi.make_track_interpolator(ra_g, dec_g, v_g, sm)

    def model_for(daz, dele):
        ra_t, dec_t = velocity.pointing_radec(tb, cfg, daz, dele)
        return ev(ra_t, dec_t, v)

    def grid_search(center, half_range, step):
        offs = np.arange(-half_range, half_range + step / 2, step)
        R = np.full((len(offs), len(offs)), -np.inf)
        best = (-np.inf, 0.0, 0.0, np.nan)
        for i, da in enumerate(offs + center[0]):
            for j, de in enumerate(offs + center[1]):
                r, g = _score(O, model_for(da, de))
                R[i, j] = r
                if r > best[0]:
                    best = (r, da, de, g)
        return best, offs + center[0], offs + center[1], R

    print("[compare] coarse pointing search ...")
    best, az1, el1, R1 = grid_search((0.0, 0.0), float(cc["search_deg"]),
                                     float(cc["coarse_step_deg"]))
    print(f"[compare]   coarse best r={best[0]:.3f} at "
          f"daz={best[1]:+.2f}, del={best[2]:+.2f}")
    print("[compare] fine search ...")
    best, az2, el2, R2 = grid_search((best[1], best[2]),
                                     float(cc["coarse_step_deg"]),
                                     float(cc["fine_step_deg"]))
    r, daz, dele, gain = best
    print(f"[compare]   fine best r={r:.3f} at daz={daz:+.2f}, del={dele:+.2f}, "
          f"gain={gain:.3f}")
    return {
        "r": r, "daz_deg": daz, "del_deg": dele, "gain": gain,
        "grids": ((az1, el1, R1), (az2, el2, R2)),
        "tb": tb, "O": O, "model_best": model_for(daz, dele),
        "model_assumed": model_for(0.0, 0.0),
    }


def make_plots(prod: pathlib.Path, cfg, v, fit):
    tb, O = fit["tb"], fit["O"]
    th = (tb - tb[0]) / 3600.0
    Mb, Ma = fit["model_best"], fit["model_assumed"]
    g = fit["gain"] if np.isfinite(fit["gain"]) and fit["gain"] > 0 else 1.0

    vm = np.nanpercentile(np.concatenate([O / g, Mb]), [2, 99.5])
    fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
    for ax, D, title in [
        (axes[0], O / g, "observed (rescaled)"),
        (axes[1], Mb, f"HI4PI @ fitted pointing\n(daz={fit['daz_deg']:+.1f}, "
                      f"del={fit['del_deg']:+.1f} deg)"),
        (axes[2], Ma, "HI4PI @ assumed pointing"),
    ]:
        im = ax.pcolormesh(v, th, D, vmin=vm[0], vmax=vm[1], cmap="inferno",
                           shading="nearest")
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("$v_{LSR}$ [km/s]")
    axes[0].set_ylabel("hours since start")
    fig.colorbar(im, ax=axes, label="$T_B$ [K]", shrink=0.85)
    fig.savefig(prod / "hi4pi_side_by_side.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # overlay: observed image + model contours
    fig, ax = plt.subplots(figsize=(9, 7))
    im = ax.pcolormesh(v, th, O / g, vmin=vm[0], vmax=vm[1], cmap="gray_r",
                       shading="nearest")
    levels = np.nanpercentile(Mb, [60, 75, 87, 95, 99])
    levels = np.unique(levels[levels > 0])
    if levels.size >= 2:
        ax.contour(v, th, np.nan_to_num(Mb), levels=levels, cmap="autumn",
                   linewidths=1.0)
    ax.set_xlabel("$v_{LSR}$ [km/s]")
    ax.set_ylabel("hours since start")
    ax.set_title("observed (grayscale) with HI4PI contours at fitted pointing")
    fig.colorbar(im, ax=ax, label="$T_B$ [K]")
    fig.tight_layout()
    fig.savefig(prod / "hi4pi_overlay.png", dpi=130)
    plt.close(fig)

    # correlation maps
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, (azg, elg, R), title in zip(
            axes, fit["grids"], ["coarse search", "fine search"]):
        im = ax.pcolormesh(elg, azg, R, cmap="viridis", shading="nearest")
        ax.plot(fit["del_deg"], fit["daz_deg"], "r+", ms=14, mew=2)
        ax.set_xlabel("elevation offset [deg]")
        ax.set_ylabel("azimuth offset [deg]")
        ax.set_title(f"{title}: correlation with HI4PI")
        fig.colorbar(im, ax=ax, label="pearson r")
    fig.tight_layout()
    fig.savefig(prod / "hi4pi_pointing_fit.png", dpi=130)
    plt.close(fig)

    # spectra overlays at a few times
    n = len(tb)
    picks = np.unique(np.linspace(0, n - 1, 4).astype(int))
    fig, axes = plt.subplots(len(picks), 1, figsize=(9, 2.6 * len(picks)),
                             sharex=True)
    axes = np.atleast_1d(axes)
    for ax, k in zip(axes, picks):
        ax.plot(v, O[k] / g, color="k", lw=1.0, label="observed (rescaled)")
        ax.plot(v, Mb[k], color="crimson", lw=1.2, alpha=0.8,
                label="HI4PI, fitted pointing")
        ax.set_ylabel("$T_B$ [K]")
        ax.text(0.02, 0.85, f"t = {th[k]:.1f} h", transform=ax.transAxes)
        ax.grid(alpha=0.3)
    axes[0].legend(loc="upper right", fontsize=8)
    axes[-1].set_xlabel("$v_{LSR}$ [km/s]")
    fig.suptitle("observed vs HI4PI spectra", y=0.995)
    fig.tight_layout()
    fig.savefig(prod / "hi4pi_spectra_overlay.png", dpi=130)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fit pointing & calibration against HI4PI")
    ap.add_argument("session", help="session directory (must be hi-processed)")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--strip", default=None,
                    help="strip cache path (default: <hi4pi.cache_dir>/strip_cache.npz)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    session_dir = pathlib.Path(args.session)
    strip = pathlib.Path(args.strip) if args.strip else \
        pathlib.Path(cfg["hi4pi"]["cache_dir"]) / "strip_cache.npz"
    if not strip.exists():
        raise FileNotFoundError(f"{strip} not found -- run hi-fetch-hi4pi first")

    t, v, T = _load_calibrated(session_dir)
    fit = fit_pointing(cfg, strip, t, v, T)

    p = cfg["pointing"]
    tsys_true = float(cfg["processing"]["tsys_assumed_k"]) / fit["gain"] \
        if np.isfinite(fit["gain"]) and fit["gain"] > 0 else None
    result = {
        "pearson_r": fit["r"],
        "daz_deg": fit["daz_deg"],
        "del_deg": fit["del_deg"],
        "az_fitted_deg": p["az_deg"] + fit["daz_deg"],
        "el_fitted_deg": p["el_deg"] + fit["del_deg"],
        "gain_obs_over_model": fit["gain"],
        "tsys_implied_k": tsys_true,
        "strip_cache": str(strip),
    }
    prod = session_dir / "products"
    prod.mkdir(exist_ok=True)
    with open(prod / "pointing_fit.json", "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
    make_plots(prod, cfg, v, fit)

    print("\n[compare] ================= RESULT =================")
    print(f"[compare] correlation with HI4PI: r = {fit['r']:.3f}")
    print(f"[compare] pointing offset: daz = {fit['daz_deg']:+.2f} deg, "
          f"del = {fit['del_deg']:+.2f} deg")
    print(f"[compare] fitted pointing: az = {result['az_fitted_deg']:.2f}, "
          f"el = {result['el_fitted_deg']:.2f} deg")
    if tsys_true:
        print(f"[compare] implied T_sys = {tsys_true:.0f} K "
              f"(update processing.tsys_assumed_k and re-run hi-process "
              f"for a true Kelvin scale)")
    print(f"[compare] plots + pointing_fit.json in {prod}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
