"""Analyze a single HI spectrum from a CSV file (e.g. from the PICTOR
web telescope, or any frequency/power table).

Runs the same science stages as the drift-scan pipeline on one spectrum:
velocity axis -> optional LSR correction -> sigma-clipped baseline ->
noise & emission mask -> moments -> multi-Gaussian decomposition -> and,
if an HI4PI cache covers the pointing, a beam-matched model overlay that
also calibrates the (arbitrary) power units to Kelvin, enabling N_HI.

Examples:
    hi-spectrum obs.csv --pictor --l 120 --b 5 --time 2026-07-04T21:30:00
    hi-spectrum obs.csv --ra 200 --dec 40 --fetch-hi4pi
CSV parsing is forgiving: comments/headers are skipped, the first numeric
column is frequency (Hz/MHz/GHz auto-detected), the last is power. A
single-column file works too if you pass --center-mhz/--bandwidth-mhz.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import calibrate, extract, hi4pi
from .config import load_config
from .velocity import F0_HI, v_topo_kms

# PICTOR is in Athens, Greece (approximate coordinates -- fine for LSR,
# where 100 km of site error changes the correction by < 0.05 km/s)
PICTOR_SITE = (38.0, 23.7, 100.0)
PICTOR_BEAM_FWHM = 10.0  # ~1.22*lambda/D for a 1.5 m dish


def read_csv_spectrum(path: pathlib.Path, freq_col: int | None,
                      power_col: int | None, center_mhz: float | None,
                      bandwidth_mhz: float | None):
    """Return (freq_hz, power) from a loosely formatted CSV/TSV."""
    rows = []
    ncol_expect = None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith(("#", "%", "//")):
            continue
        for delim in (",", ";", "\t", None):
            parts = line.split(delim)
            try:
                vals = [float(p) for p in parts if p != ""]
            except ValueError:
                continue
            if vals:
                if ncol_expect is None:
                    ncol_expect = len(vals)
                if len(vals) == ncol_expect:
                    rows.append(vals)
                break
    if not rows:
        raise SystemExit(f"[spectrum] no numeric rows found in {path}")
    data = np.array(rows)
    print(f"[spectrum] {path.name}: {data.shape[0]} rows x {data.shape[1]} columns")

    if data.shape[1] == 1:
        if center_mhz is None or bandwidth_mhz is None:
            raise SystemExit("[spectrum] single-column file: pass --center-mhz "
                             "and --bandwidth-mhz to define the frequency axis")
        n = data.shape[0]
        freq_hz = (center_mhz + bandwidth_mhz * (np.arange(n) / (n - 1) - 0.5)) * 1e6
        return freq_hz, data[:, 0]

    fc = 0 if freq_col is None else freq_col
    pc = data.shape[1] - 1 if power_col is None else power_col
    f = data[:, fc]
    med = np.median(np.abs(f))
    if med > 1e8:
        freq_hz = f                       # Hz
    elif med > 1e2:
        freq_hz = f * 1e6                 # MHz
    else:
        freq_hz = f * 1e9                 # GHz
    if not (1.3e9 < np.median(freq_hz) < 1.6e9):
        print("[spectrum] WARNING: frequency axis is far from 1420 MHz -- "
              "check units / --freq-col")
    return freq_hz, data[:, pc]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Analyze a single HI spectrum CSV")
    ap.add_argument("csv", help="spectrum file (e.g. from PICTOR)")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--out", default=None, help="output dir (default <csv>_analysis)")
    ap.add_argument("--freq-col", type=int, default=None)
    ap.add_argument("--power-col", type=int, default=None)
    ap.add_argument("--center-mhz", type=float, default=None)
    ap.add_argument("--bandwidth-mhz", type=float, default=None)
    # pointing (either equatorial or galactic)
    ap.add_argument("--ra", type=float, default=None)
    ap.add_argument("--dec", type=float, default=None)
    ap.add_argument("--l", type=float, default=None, dest="gal_l")
    ap.add_argument("--b", type=float, default=None, dest="gal_b")
    # observation circumstances (for the LSR correction)
    ap.add_argument("--time", default=None, help="UTC, e.g. 2026-07-04T21:30:00")
    ap.add_argument("--site-lat", type=float, default=None)
    ap.add_argument("--site-lon", type=float, default=None)
    ap.add_argument("--pictor", action="store_true",
                    help="use the PICTOR (Athens) site and a 10 deg beam")
    ap.add_argument("--beam-fwhm", type=float, default=None)
    # HI4PI comparison
    ap.add_argument("--strip", default=None, help="existing strip cache .npz")
    ap.add_argument("--fetch-hi4pi", action="store_true",
                    help="download the 1-4 HI4PI tiles around the pointing")
    ap.add_argument("--baseline-order", type=int, default=3)
    ap.add_argument("--edge-frac", type=float, default=0.05,
                    help="fraction of channels dropped at each band edge "
                         "before baseline fitting (filter rolloff region)")
    args = ap.parse_args(argv)

    csv_path = pathlib.Path(args.csv)
    out = pathlib.Path(args.out) if args.out else \
        csv_path.with_name(csv_path.stem + "_analysis")
    out.mkdir(parents=True, exist_ok=True)

    # ---- pointing ---------------------------------------------------------------
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    sc = None
    if args.gal_l is not None and args.gal_b is not None:
        sc = SkyCoord(l=args.gal_l * u.deg, b=args.gal_b * u.deg,
                      frame="galactic").icrs
    elif args.ra is not None and args.dec is not None:
        sc = SkyCoord(args.ra * u.deg, args.dec * u.deg)
    if sc is not None:
        print(f"[spectrum] pointing: RA {sc.ra.deg:.2f}, Dec {sc.dec.deg:+.2f} "
              f"(l={sc.galactic.l.deg:.1f}, b={sc.galactic.b.deg:+.1f})")

    # ---- read + velocity axis -----------------------------------------------------
    freq_hz, power = read_csv_spectrum(csv_path, args.freq_col, args.power_col,
                                       args.center_mhz, args.bandwidth_mhz)
    if freq_hz[0] > freq_hz[-1]:
        freq_hz, power = freq_hz[::-1], power[::-1]
    ne = int(args.edge_frac * len(freq_hz))
    if ne > 0:
        freq_hz, power = freq_hz[ne:-ne], power[ne:-ne]
    v = v_topo_kms(freq_hz)  # descending when freq ascending

    vcorr = 0.0
    site = None
    if args.pictor:
        site = PICTOR_SITE
    if args.site_lat is not None and args.site_lon is not None:
        site = (args.site_lat, args.site_lon, 100.0)
    if args.time and sc is not None and site is not None:
        from astropy.coordinates import EarthLocation
        from astropy.time import Time
        from .velocity import _LSR_APEX, _LSR_SPEED
        loc = EarthLocation(lat=site[0] * u.deg, lon=site[1] * u.deg,
                            height=site[2] * u.m)
        tt = Time(args.time, scale="utc")
        vb = sc.radial_velocity_correction(kind="barycentric", obstime=tt,
                                           location=loc).to_value(u.km / u.s)
        vcorr = float(vb + _LSR_SPEED * np.cos(sc.separation(_LSR_APEX).rad))
        print(f"[spectrum] LSR correction: {vcorr:+.2f} km/s")
    else:
        print("[spectrum] no --time/--site/pointing given -> velocities are "
              "TOPOCENTRIC (can be off by up to ~30 km/s from v_LSR)")
    v_lsr = v + vcorr

    # ---- baseline + stats (in ascending-velocity order) ----------------------------
    order = np.argsort(v_lsr)
    vv, ss = v_lsr[order], power[order].astype(float)
    base, _ = calibrate.fit_baseline(ss, vv, order=args.baseline_order)
    line = ss - base
    sig = extract.noise_sigma(line, vv)
    mask = extract.emission_mask(line, sig)
    comps = extract.decompose(vv, line, sig)
    print(f"[spectrum] noise sigma = {sig:.3g} (input units), "
          f"{int(mask.sum())} emission channels, {len(comps)} gaussian components")
    for i, c in enumerate(comps, 1):
        print(f"[spectrum]   G{i}: amp={c['amp_k']:.3g}  v0={c['v0_kms']:+7.1f} km/s"
              f"  sigma_v={c['sigma_v_kms']:5.1f} km/s")

    # ---- optional HI4PI model + Kelvin calibration ---------------------------------
    model = None
    kelvin_per_unit = None
    beam = args.beam_fwhm or (PICTOR_BEAM_FWHM if args.pictor
                              else load_config(args.config)["pointing"]["beam_fwhm_deg"])
    strip_path = pathlib.Path(args.strip) if args.strip else None
    if args.fetch_hi4pi and sc is not None:
        cfg = load_config(args.config)
        strip_path = pathlib.Path(cfg["hi4pi"]["cache_dir"]) / \
            f"point_{sc.ra.deg:.0f}_{sc.dec.deg:+.0f}.npz"
        if not strip_path.exists():
            m = beam * 1.2
            hi4pi.build_strip_cache(cfg, sc.dec.deg - m, sc.dec.deg + m,
                                    out_path=strip_path,
                                    ra_min=(sc.ra.deg - m) % 360,
                                    ra_max=(sc.ra.deg + m) % 360)
    if strip_path is not None and strip_path.exists() and sc is not None:
        ra_g, dec_g, v_g, cube = hi4pi.load_strip(strip_path)
        sm = hi4pi.beam_smooth_strip(ra_g, dec_g, cube, beam)
        model = hi4pi.track_spectra(ra_g, dec_g, v_g, sm,
                                    [sc.ra.deg], [sc.dec.deg], vv)[0]
        fin = np.isfinite(model) & mask
        if fin.sum() > 5 and np.nansum(line[fin] ** 2) > 0:
            kelvin_per_unit = float(np.nansum(line[fin] * model[fin])
                                    / np.nansum(line[fin] ** 2))
            dv = float(np.median(np.diff(vv)))
            mom0_k = float(np.sum(line[mask]) * dv * kelvin_per_unit)
            print(f"[spectrum] HI4PI amplitude fit: 1 input unit = "
                  f"{kelvin_per_unit:.3g} K")
            print(f"[spectrum] integrated line: {mom0_k:.1f} K km/s -> "
                  f"N_HI = {1.823e18 * mom0_k:.3g} cm^-2")
        else:
            print("[spectrum] pointing outside cache or no emission overlap -- "
                  "no Kelvin calibration")

    # ---- outputs -------------------------------------------------------------------
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(9, 8))
    a1.plot(freq_hz / 1e6, power, lw=0.8, color="k")
    a1.plot(freq_hz / 1e6, base[np.argsort(order)], lw=1.0, color="orange",
            label="fitted baseline")
    a1.axvline(F0_HI / 1e6, color="r", ls="--", lw=0.8, label="HI rest")
    a1.set_xlabel("frequency [MHz]")
    a1.set_ylabel("power [input units]")
    a1.set_title(csv_path.name)
    a1.legend(fontsize=8)
    a1.grid(alpha=0.3)

    scale = kelvin_per_unit or 1.0
    a2.plot(vv, line * scale, lw=0.9, color="k", label="data - baseline")
    if comps:
        fit = extract._multi_gauss(vv, *[p for c in comps for p in
                                         (c["amp_k"], c["v0_kms"], c["sigma_v_kms"])])
        a2.plot(vv, fit * scale, lw=1.2, color="royalblue", alpha=0.8,
                label=f"{len(comps)}-gaussian fit")
    if model is not None:
        a2.plot(vv, model, lw=1.2, color="crimson", alpha=0.8,
                label=f"HI4PI x {beam:.0f} deg beam")
    a2.set_xlabel("$v_{LSR}$ [km/s]" if vcorr or model is not None
                  else "$v_{topo}$ [km/s]")
    a2.set_ylabel("$T_B$ [K]" if kelvin_per_unit else "line [input units]")
    a2.legend(fontsize=8)
    a2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "spectrum.png", dpi=130)
    plt.close(fig)

    summary = {
        "file": str(csv_path), "n_channels": int(len(vv)),
        "vcorr_lsr_kms": vcorr, "noise_sigma_input_units": sig,
        "kelvin_per_unit": kelvin_per_unit, "beam_fwhm_deg": beam,
        "components": comps,
    }
    with open(out / "summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2)
    print(f"[spectrum] wrote {out / 'spectrum.png'} and summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
