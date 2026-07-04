"""End-to-end processing of a raw session:

raw ON/OFF chunks -> bandpass calibration -> RFI flagging -> LSR regridding
-> baseline removal -> time & RA waterfalls -> science extraction -> plots.

Usage:  hi-process data/raw/20260704_120000_obs [--config config.yaml]
Products land in <session>/products/.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import calibrate, extract, rfi, session, velocity
from .config import load_config


def process_session(session_dir: pathlib.Path, cfg: dict,
                    start_h: float | None = None, stop_h: float | None = None,
                    last_h: float | None = None) -> pathlib.Path:
    pc = cfg["processing"]
    meta, t, on, off = session.load_session(session_dir)
    nfft = int(meta["nfft"])

    # optional time-window selection (hours relative to session start)
    rel_h = (t - t[0]) / 3600.0
    sel = np.ones(len(t), dtype=bool)
    if last_h is not None:
        sel &= rel_h >= rel_h[-1] - last_h
    if start_h is not None:
        sel &= rel_h >= start_h
    if stop_h is not None:
        sel &= rel_h <= stop_h
    if not sel.all():
        print(f"[process] time window keeps {sel.sum()}/{len(t)} cycles")
        t, on, off = t[sel], on[sel], off[sel]
    if len(t) < 4:
        raise SystemExit("[process] fewer than 4 cycles selected -- nothing to do")
    print(f"[process] {session_dir}: {len(t)} cycles, "
          f"{(t[-1]-t[0])/3600:.2f} h span")

    prod = session_dir / "products"
    prod.mkdir(exist_ok=True)

    # --- calibration to Kelvin ------------------------------------------------
    tsys = float(pc["tsys_assumed_k"])
    q = calibrate.quotient_kelvin(on, off, tsys)
    good_chan = calibrate.channel_mask(nfft, float(pc["edge_frac"]),
                                       int(pc["dc_halfwidth_bins"]))

    freqs = velocity.channel_freqs_hz(float(meta["freq_on_hz"]),
                                      float(meta["sample_rate_hz"]), nfft)
    v_topo = velocity.v_topo_kms(freqs)
    protect = np.abs(v_topo) < 220.0  # galactic HI can live here

    bad, bad_cycles = rfi.flag_rfi(
        q, good_chan, protect,
        cycle_sigma=float(pc["rfi_cycle_sigma"]),
        chan_sigma=float(pc["rfi_chan_sigma"]),
        pixel_sigma=float(pc["rfi_pixel_sigma"]))
    q[bad] = np.nan
    keep = ~bad_cycles
    q, t = q[keep], t[keep]

    # --- pointing track and LSR velocities -------------------------------------
    ra, dec = velocity.pointing_radec(t, cfg)
    vcorr = velocity.vlsr_correction_kms(t, ra, dec, cfg)
    print(f"[process] track: dec ~ {np.mean(dec):+.2f} deg, "
          f"RA {ra.min():.1f}..{ra.max():.1f} deg, "
          f"v_LSR corr {vcorr.min():+.1f}..{vcorr.max():+.1f} km/s")

    vgrid = np.arange(float(pc["v_min_kms"]),
                      float(pc["v_max_kms"]) + 1e-6, float(pc["dv_kms"]))
    # v_topo decreases with frequency; flip to ascending for np.interp
    T = np.full((len(t), len(vgrid)), np.nan, dtype=np.float32)
    for i in range(len(t)):
        fin = np.isfinite(q[i])
        if fin.sum() < 32:
            continue
        x = (v_topo + vcorr[i])[fin][::-1]
        y = q[i][fin][::-1]
        T[i] = np.interp(vgrid, x, y, left=np.nan, right=np.nan)

    # --- per-spectrum polynomial baselines --------------------------------------
    T = calibrate.remove_baselines(T, vgrid, order=int(pc["baseline_order"]))

    np.savez_compressed(prod / "calibrated.npz",
                        t=t, v=vgrid, T=T, ra=ra, dec=dec, vcorr=vcorr,
                        meta=json.dumps(meta))
    print(f"[process] calibrated waterfall saved: {prod/'calibrated.npz'}")

    # --- RA-binned map ----------------------------------------------------------
    ra_step = float(pc["ra_bin_deg"])
    edges = np.arange(0.0, 360.0 + ra_step, ra_step)
    idx = np.digitize(ra, edges) - 1
    nb = len(edges) - 1
    W = np.full((nb, len(vgrid)), np.nan, dtype=np.float32)
    dec_b = np.full(nb, np.nan)
    counts = np.zeros(nb, dtype=int)
    for bidx in range(nb):
        sel = idx == bidx
        if sel.sum() == 0:
            continue
        W[bidx] = np.nanmean(T[sel], axis=0)
        dec_b[bidx] = np.mean(dec[sel])
        counts[bidx] = sel.sum()
    filled = counts > 0
    ra_b = 0.5 * (edges[:-1] + edges[1:])
    l_b = np.full(nb, np.nan)
    b_b = np.full(nb, np.nan)
    if filled.any():
        l_b[filled], b_b[filled] = velocity.galactic_lb(ra_b[filled], dec_b[filled])
    np.savez_compressed(prod / "ra_map.npz", ra=ra_b, dec=dec_b, l=l_b, b=b_b,
                        v=vgrid, W=W, counts=counts)

    _plots(prod, t, vgrid, T, ra_b, W, filled, meta)

    # --- science extraction -------------------------------------------------------
    rows = extract.analyze_bins(vgrid, W[filled], ra_b[filled], dec_b[filled],
                                l_b[filled], b_b[filled], prod / "extraction.csv")
    _extraction_plots(prod, vgrid, W[filled], ra_b[filled], l_b[filled],
                      b_b[filled], rows)
    return prod


def _plots(prod, t, v, T, ra_b, W, filled, meta):
    th = (t - t[0]) / 3600.0
    fig, ax = plt.subplots(figsize=(10, 6))
    vm = np.nanpercentile(T, [2, 99.5])
    im = ax.pcolormesh(v, th, T, vmin=vm[0], vmax=vm[1], cmap="inferno",
                       shading="nearest")
    ax.set_xlabel("$v_{LSR}$ [km/s]")
    ax.set_ylabel("hours since start")
    ax.set_title(f"drift-scan waterfall -- {meta.get('kind','obs')}")
    fig.colorbar(im, ax=ax, label="$T_B$ [K]")
    fig.tight_layout()
    fig.savefig(prod / "waterfall_time.png", dpi=130)
    plt.close(fig)

    if filled.sum() > 1:
        fig, ax = plt.subplots(figsize=(11, 5))
        Wm = np.where(np.isfinite(W), W, np.nan)
        im = ax.pcolormesh(ra_b, v, Wm.T, vmin=vm[0], vmax=vm[1],
                           cmap="inferno", shading="nearest")
        ax.set_xlabel("RA [deg]")
        ax.set_ylabel("$v_{LSR}$ [km/s]")
        ax.set_title("RA x velocity map (this is your slice of the Milky Way)")
        fig.colorbar(im, ax=ax, label="$T_B$ [K]")
        fig.tight_layout()
        fig.savefig(prod / "waterfall_ra.png", dpi=130)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    avg = np.nanmean(T, axis=0)
    ax.plot(v, avg, lw=1.0, color="k")
    ax.axvline(0, color="0.7", lw=0.7)
    ax.set_xlabel("$v_{LSR}$ [km/s]")
    ax.set_ylabel("$T_B$ [K]")
    ax.set_title("session-averaged spectrum")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(prod / "avg_spectrum.png", dpi=130)
    plt.close(fig)


def _extraction_plots(prod, v, W, ra_b, l_b, b_b, rows):
    if rows:
        ra = np.array([r["ra_deg"] for r in rows])
        nhi = np.array([r["nhi_cm2"] for r in rows])
        pk = np.array([r["peak_k"] for r in rows])
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        a1.plot(ra, nhi, ".-", ms=3)
        a1.set_ylabel(r"$N_{HI}$ [cm$^{-2}$]")
        a1.grid(alpha=0.3)
        a2.plot(ra, pk, ".-", ms=3, color="crimson")
        a2.set_ylabel("peak $T_B$ [K]")
        a2.set_xlabel("RA [deg]")
        a2.grid(alpha=0.3)
        a1.set_title("HI column density and peak brightness along the drift strip")
        fig.tight_layout()
        fig.savefig(prod / "column_density.png", dpi=130)
        plt.close(fig)

    xs, ys, ts = extract.faceon_points(v, W, l_b, b_b)
    if xs.size > 20:
        fig, ax = plt.subplots(figsize=(7.5, 7.5))
        sc = ax.scatter(xs, ys, c=ts, s=6, cmap="viridis",
                        vmin=0, vmax=np.percentile(ts, 98))
        ax.plot(0, extract.R0_KPC, "*", color="orange", ms=14, label="Sun")
        ax.plot(0, 0, "k+", ms=12, label="Gal. center")
        ax.set_xlabel("x [kpc]")
        ax.set_ylabel("y [kpc]")
        ax.set_title("face-on HI map from kinematic distances (outer Galaxy)")
        ax.set_aspect("equal")
        ax.legend(loc="lower right")
        fig.colorbar(sc, ax=ax, label="$T_B$ [K]")
        fig.tight_layout()
        fig.savefig(prod / "faceon_map.png", dpi=130)
        plt.close(fig)
        print(f"[extract] face-on map: {xs.size} kinematic points")

    rc = extract.rotation_curve(v, W, l_b, b_b)
    if rc.shape[0] > 3:
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(rc[:, 0], rc[:, 1], "o", ms=4)
        ax.set_xlabel("R [kpc]")
        ax.set_ylabel("$V_{c}$ [km/s]")
        ax.set_title("tangent-point rotation curve (inner Galaxy)")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(prod / "rotation_curve.png", dpi=130)
        plt.close(fig)
        np.savetxt(prod / "rotation_curve.csv", rc, delimiter=",",
                   header="R_kpc,Vc_kms,l_deg", comments="")
        print(f"[extract] rotation curve: {rc.shape[0]} tangent points")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Process a raw drift-scan session")
    ap.add_argument("session", help="session directory (data/raw/...)")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--start-h", type=float, default=None,
                    help="process only cycles after this many hours into the session")
    ap.add_argument("--stop-h", type=float, default=None,
                    help="process only cycles before this many hours into the session")
    ap.add_argument("--last-h", type=float, default=None,
                    help="process only the most recent N hours (quick look)")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)
    prod = process_session(pathlib.Path(args.session), cfg,
                           start_h=args.start_h, stop_h=args.stop_h,
                           last_h=args.last_h)
    print(f"[process] all products in {prod}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
