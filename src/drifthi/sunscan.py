"""Solar pointing calibration: find where the dish actually points by
letting the Sun drift through the beam.

The Sun is a ~70 K antenna-temperature signal for an 80 cm dish -- it
roughly doubles the broadband power, so a normal `hi-observe` session that
spans the daytime hours contains a huge bump in total power. This tool:

  1. fits (baseline + Gaussian) to total power vs time,
  2. converts the peak time to the Sun's az/el (astropy),
  3. reports the pointing constraint: at closest approach the beam center
     lies on the line through the Sun PERPENDICULAR to the Sun's path, so
     the along-path component of your pointing error is measured exactly;
     the cross-path component is not (a miss above or below the path only
     lowers the bump amplitude),
  4. measures the beam FWHM from the bump duration (a free bonus: Gaussian
     beam x moving point source = Gaussian in time with the same width).

Usage (leave the dish exactly as it will be for night observations!):
    hi-observe --tag sunscan       # run from ~3 h before to ~3 h after
                                   # the Sun's closest approach (~solar noon)
    hi-sunscan data/raw/<session>_sunscan

Note: works only when the Sun's declination is within ~a beam of the dish
declination (in summer, dec_sun ~ +23: fine for el 50-60 from Lithuania).
"""

from __future__ import annotations

import argparse
import json
import pathlib

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from . import session
from .calibrate import channel_mask
from .config import load_config


def _sun_altaz(t_unix, cfg):
    import astropy.units as u
    from astropy.coordinates import AltAz, get_sun
    from astropy.time import Time
    from .velocity import earth_location

    t = Time(np.atleast_1d(t_unix), format="unix")
    aa = get_sun(t).transform_to(AltAz(obstime=t, location=earth_location(cfg)))
    return aa.az.deg, aa.alt.deg


def _model(t, base, slope, amp, t0, sig):
    return base + slope * (t - t0) + amp * np.exp(-0.5 * ((t - t0) / sig) ** 2)


def _predict(cfg) -> int:
    """When does the Sun cross the configured pointing, and how strongly?"""
    import time as _time
    pnt = cfg["pointing"]
    fwhm = float(pnt["beam_fwhm_deg"])
    ts = _time.time() + np.arange(0, 48 * 3600, 120.0)
    az, el = _sun_altaz(ts, cfg)
    ce = np.cos(np.radians(el))
    dx = (az - float(pnt["az_deg"]) + 180.0) % 360.0 - 180.0
    sep = np.hypot(dx * ce, el - float(pnt["el_deg"]))
    up = el > 0
    if not up.any():
        print("[sunscan] the Sun never rises in the next 48 h?!")
        return 1
    for day, sl in (("next 24 h", slice(0, 720)), ("following day", slice(720, 1440))):
        s = sep[sl].copy()
        s[~up[sl]] = np.inf
        i = int(np.argmin(s))
        if not np.isfinite(s[i]):
            continue
        t_best = ts[sl][i]
        amp = np.exp(-4 * np.log(2) * (s[i] / fwhm) ** 2)
        verdict = ("STRONG transit -- go" if amp > 0.5 else
                   "detectable" if amp > 0.05 else
                   "too far from the beam -- no usable transit")
        print(f"[sunscan] {day}: closest approach "
              f"{_time.strftime('%Y-%m-%d %H:%M', _time.gmtime(t_best))} UTC, "
              f"separation {s[i]:.1f} deg, bump ~{100*amp:.0f}% of max -- {verdict}")
    print(f"[sunscan] record from ~2.5 h before to ~2.5 h after the closest approach")
    return 0


def _sun_at(cfg, when_utc: str) -> int:
    """Where will the Sun be at a given time? Point the dish there and the
    transit happens exactly then, dead-center."""
    from astropy.time import Time
    t = Time(when_utc, scale="utc").unix
    az, el = (float(x[0]) for x in _sun_altaz(np.array([t]), cfg))
    if el < 5:
        print(f"[sunscan] Sun is at el={el:.1f} deg then -- below/near the "
              "horizon, pick another time")
        return 1
    from .velocity import pointing_radec
    cfg2 = {**cfg, "pointing": {**cfg["pointing"], "az_deg": az, "el_deg": el}}
    _, dec = pointing_radec(np.array([t]), cfg2)
    print(f"[sunscan] at {when_utc} UTC the Sun is at az={az:.2f}, el={el:.2f}")
    print(f"[sunscan] point the dish there (tip: align by the feed shadow), set "
          f"pointing.az_deg/el_deg in config.yaml, and you get a dead-center "
          f"transit at that time")
    print(f"[sunscan] kept for the night, that pointing drift-scans the "
          f"dec = {float(dec[0]):+.1f} deg ring")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fit a solar transit to calibrate pointing")
    ap.add_argument("session", nargs="?", default=None,
                    help="daytime session directory (omit with --predict/--sun-at)")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--predict", action="store_true",
                    help="just report when the Sun crosses the configured "
                         "pointing in the next 48 h (no data needed)")
    ap.add_argument("--sun-at", default=None, metavar="UTC",
                    help="where is the Sun at this UTC time (e.g. "
                         "2026-07-06T14:00:00)? Point the dish there to get a "
                         "dead-center transit exactly then")
    args = ap.parse_args(argv)
    cfg = load_config(args.config)

    if args.sun_at:
        return _sun_at(cfg, args.sun_at)
    if args.predict:
        return _predict(cfg)
    if args.session is None:
        ap.error("give a session directory, or use --predict / --sun-at")

    meta, t, on, off = session.load_session(pathlib.Path(args.session))
    good = channel_mask(int(meta["nfft"]), 0.10, 3)
    p = 0.5 * (np.nanmedian(on[:, good], axis=1) + np.nanmedian(off[:, good], axis=1))
    pn = p / np.median(p)

    from scipy.ndimage import median_filter
    from scipy.optimize import curve_fit
    sm = median_filter(pn, size=max(3, len(pn) // 100 | 1))
    i0 = int(np.argmax(sm))
    p0 = [float(np.median(pn)), 0.0, max(float(sm[i0] - np.median(pn)), 0.01),
          float(t[i0]), 1500.0]
    span = t[-1] - t[0]
    bounds = ([0.0, -1e-5, 0.0, t[0], 120.0], [10.0, 1e-5, 10.0, t[-1], span])
    popt, _ = curve_fit(_model, t, pn, p0=p0, bounds=bounds, maxfev=20000)
    base, slope, amp, t0, sig = popt
    resid = pn - _model(t, *popt)
    noise = 1.4826 * np.median(np.abs(resid - np.median(resid)))
    snr = amp / max(noise, 1e-9)
    print(f"[sunscan] transit fit: amplitude {100*amp:.1f}% of system power, "
          f"SNR {snr:.0f}, sigma_t {sig:.0f} s")
    if snr < 5:
        print("[sunscan] WARNING: no significant solar bump -- the Sun may not "
              "cross your beam (dish dec vs sun dec), or the session missed "
              "the transit time")

    # Sun geometry at closest approach, and its path direction
    az_s, el_s = (float(x[0]) for x in _sun_altaz(np.array([t0]), cfg))
    az_p, el_p = _sun_altaz(np.array([t0 - 60, t0 + 60]), cfg)
    # tangent-plane vectors: x = az*cos(el), y = el
    ce = np.cos(np.radians(el_s))
    u = np.array([(az_p[1] - az_p[0]) * ce, el_p[1] - el_p[0]])
    omega = np.linalg.norm(u) / 120.0            # deg/s on the sky
    u /= np.linalg.norm(u)
    n = np.array([-u[1], u[0]])                  # perpendicular to the path
    fwhm = 2.3548 * sig * omega

    pnt = cfg["pointing"]
    az0, el0 = float(pnt["az_deg"]), float(pnt["el_deg"])
    e = np.array([(az0 - az_s) * ce, el0 - el_s])   # assumed - sun, on sky
    along = float(e @ u)                             # measured error component
    az_corr = az0 - along * u[0] / ce
    el_corr = el0 - along * u[1]

    import time as _time
    iso = _time.strftime("%Y-%m-%d %H:%M:%S", _time.gmtime(t0))
    print(f"[sunscan] closest approach {iso} UTC; Sun at az={az_s:.2f}, el={el_s:.2f}")
    print(f"[sunscan] measured beam FWHM ~ {fwhm:.1f} deg "
          f"(config says {pnt['beam_fwhm_deg']}) -- update config if trustworthy")
    print(f"[sunscan] along-path pointing error: {along:+.2f} deg")
    print(f"[sunscan] corrected pointing (along-path component fixed): "
          f"az={az_corr:.2f}, el={el_corr:.2f}")
    print(f"[sunscan] the cross-path component is NOT constrained by one transit "
          f"(it only lowers the bump amplitude); hi-compare against HI4PI nails it, "
          f"or repeat the sunscan after tilting the dish ~5 deg in elevation")

    prod = pathlib.Path(args.session) / "products"
    prod.mkdir(exist_ok=True)
    th = (t - t[0]) / 3600.0
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(th, pn, ".", ms=2, color="0.6", label="total power (per cycle)")
    tt = np.linspace(t[0], t[-1], 800)
    ax.plot((tt - t[0]) / 3600, _model(tt, *popt), "crimson", lw=1.5, label="fit")
    ax.axvline((t0 - t[0]) / 3600, color="k", ls="--", lw=0.8)
    ax.set_xlabel("hours since session start")
    ax.set_ylabel("relative broadband power")
    ax.set_title(f"solar transit: closest approach {iso} UTC, "
                 f"beam FWHM ~ {fwhm:.1f} deg")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(prod / "sunscan.png", dpi=130)
    plt.close(fig)

    with open(prod / "sunscan.json", "w", encoding="utf-8") as fh:
        json.dump({"t0_unix": float(t0), "t0_utc": iso, "sun_az_deg": az_s,
                   "sun_el_deg": el_s, "amp_frac": float(amp), "snr": float(snr),
                   "beam_fwhm_deg": float(fwhm), "along_path_error_deg": along,
                   "az_corrected_deg": float(az_corr),
                   "el_corrected_deg": float(el_corr)}, fh, indent=2)
    print(f"[sunscan] wrote {prod/'sunscan.png'} and sunscan.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
