"""Auto-scheduler for the PICTOR web telescope (pictortelescope.com).

PICTOR is a 1.5 m dish fixed at the zenith over Athens. The galactic plane
crosses its beam twice per sidereal day: once in Cygnus (l ~ 75 deg, the
brightest HI in the northern sky, ~03:00 EEST at the moment) and once near
the anticenter in Perseus/Auriga (l ~ 165 deg). This bot computes those
crossing times with astropy, sleeps until just before one, and then submits
a chain of observation requests (the site caps each at 600 s) covering the
passage -- with raw_data=1 so every CSV lands in your inbox while you sleep.

Typical use (leave running in a terminal in the evening):
    hi-pictor-bot --email you@example.com                # next Cygnus pass
    hi-pictor-bot --email you@example.com --target anticenter --nights 3
    hi-pictor-bot --email you@example.com --dry-run      # just print the plan

Be a good guest: the default 1 h window = 6 requests per night on a shared
public instrument. Analyze each CSV with:
    hi-spectrum obs.csv --pictor --ra <ra> --dec 38.0 --time <utc>
(the bot prints the exact command, with the zenith RA, for every slot).
"""

from __future__ import annotations

import argparse
import sys
import time
import urllib.parse
import urllib.request

import numpy as np

# line-buffer stdout even when redirected to a log file (nohup ... >> log),
# so the schedule is visible immediately instead of sitting in a buffer
try:
    sys.stdout.reconfigure(line_buffering=True)
except AttributeError:
    pass

OBSERVE_URL = "https://pictortelescope.com/observe"
PICTOR_LAT, PICTOR_LON = 38.0, 23.7
MAX_DUR_S = 600


def _zenith(t_unix: np.ndarray):
    """ICRS SkyCoord of the PICTOR zenith at unix times."""
    import astropy.units as u
    from astropy.coordinates import AltAz, EarthLocation, SkyCoord
    from astropy.time import Time
    from astropy.utils import iers

    # standalone-friendly: never try to download IERS tables (sub-arcsecond
    # corrections are irrelevant for a 10 deg beam)
    iers.conf.auto_download = False
    iers.conf.auto_max_age = None

    loc = EarthLocation(lat=PICTOR_LAT * u.deg, lon=PICTOR_LON * u.deg, height=100 * u.m)
    t = Time(np.atleast_1d(t_unix), format="unix")
    return SkyCoord(AltAz(az=0 * u.deg, alt=90 * u.deg, obstime=t, location=loc)).icrs


def find_crossings(t_start_unix: float, days: float = 2.0):
    """(t_unix, l_deg, ra_deg) of galactic-plane crossings (b=0) of the zenith."""
    ts = t_start_unix + np.arange(0, days * 86400, 240.0)
    g = _zenith(ts).galactic
    b = g.b.deg
    out = []
    for i in np.flatnonzero(np.sign(b[:-1]) != np.sign(b[1:])):
        lo, hi = ts[i], ts[i + 1]
        for _ in range(20):  # bisection to ~seconds
            mid = 0.5 * (lo + hi)
            bm = _zenith(np.array([mid])).galactic.b.deg[0]
            if np.sign(bm) == np.sign(_zenith(np.array([lo])).galactic.b.deg[0]):
                lo = mid
            else:
                hi = mid
        tc = 0.5 * (lo + hi)
        zc = _zenith(np.array([tc]))
        out.append((tc, float(zc.galactic.l.deg[0]), float(zc.ra.deg[0])))
    return out


def submit_observation(name: str, email: str, duration_s: int, nbins: int,
                       f_center_mhz: float = 1420.0, dry_run: bool = False) -> bool:
    fields = {
        "obs_name": name,
        "f_center": f"{f_center_mhz:g}",
        "bandwidth": "2.4mhz",
        "channels": "2048",
        "nbins": str(nbins),
        "duration": str(int(duration_s)),
        "raw_data": "1",          # emails the CSV
        "email": email,
        "submit_btn": "",
    }
    if dry_run:
        print(f"[pictor-bot]   DRY RUN -- would POST: {fields}")
        return True
    data = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(
        OBSERVE_URL, data=data,
        headers={"User-Agent": "Mozilla/5.0 (drifthi pictor-bot)",
                 "Content-Type": "application/x-www-form-urlencoded"})
    for attempt in (1, 2):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                body = resp.read(20000).decode("utf-8", errors="replace").lower()
            ok = resp.status == 200 and ("success" in body or "queue" in body
                                         or "received" in body or "thank" in body)
            print(f"[pictor-bot]   POST -> HTTP {resp.status} "
                  f"({'looks OK' if ok else 'check your email to confirm'})")
            return True
        except OSError as e:
            print(f"[pictor-bot]   submission attempt {attempt} failed: {e!r}")
            time.sleep(10)
    return False


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Auto-observe galactic-plane transits with PICTOR")
    ap.add_argument("--email", required=True, help="where PICTOR sends the CSVs")
    ap.add_argument("--target", choices=["cygnus", "anticenter", "both"], default="cygnus",
                    help="which plane crossing to observe (cygnus = l~75, brightest)")
    ap.add_argument("--window-h", type=float, default=1.0,
                    help="observing window centered on the crossing [hours]")
    ap.add_argument("--gap-s", type=float, default=90.0,
                    help="pause between consecutive requests (queue processing time)")
    ap.add_argument("--nights", type=int, default=1, help="repeat for N nights")
    ap.add_argument("--name", default="mw_transit", help="observation name prefix")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the schedule and the would-be POSTs, submit nothing")
    args = ap.parse_args(argv)

    def is_wanted(l_deg):
        if args.target == "both":
            return True
        return (30 < l_deg < 120) == (args.target == "cygnus")

    n_per_window = max(1, int(round(args.window_h * 3600 / (MAX_DUR_S + args.gap_s))))
    print(f"[pictor-bot] plan: {args.window_h:.1f} h window -> {n_per_window} x "
          f"{MAX_DUR_S} s observations per crossing, CSVs to {args.email}")

    crossings = [c for c in find_crossings(time.time(), days=args.nights + 1.5)
                 if is_wanted(c[1])][: args.nights]
    if not crossings:
        print("[pictor-bot] no matching crossings found -- try --target both")
        return 1

    for tc, l_deg, ra_deg in crossings:
        t0 = tc - args.window_h * 1800  # window start (half-window before b=0)
        loc_str = "Cygnus" if 30 < l_deg < 120 else "Perseus/anticenter"
        print(f"\n[pictor-bot] crossing: {time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(tc))} UTC"
              f"  ({loc_str}, l={l_deg:.1f}, zenith RA={ra_deg:.1f})")
        print(f"[pictor-bot] window starts {time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(t0))} UTC")

        for k in range(n_per_window):
            t_sub = t0 + k * (MAX_DUR_S + args.gap_s)
            wait = t_sub - time.time()
            if wait > 0 and not args.dry_run:
                print(f"[pictor-bot] sleeping {wait/3600:.2f} h until slot {k+1}/{n_per_window} ...")
                time.sleep(wait)
            elif wait < -MAX_DUR_S:
                print(f"[pictor-bot] slot {k+1} already passed, skipping")
                continue
            zen = _zenith(np.array([t_sub + MAX_DUR_S / 2]))
            ra_mid = float(zen.ra.deg[0])
            stamp = time.strftime("%m%d_%H%M", time.gmtime(t_sub))
            name = f"{args.name}_{stamp}_slot{k+1}"
            print(f"[pictor-bot] slot {k+1}/{n_per_window} at "
                  f"{time.strftime('%H:%M:%S', time.gmtime(t_sub))} UTC -> submitting '{name}'")
            submit_observation(name, args.email, MAX_DUR_S, nbins=100, dry_run=args.dry_run)
            t_mid_iso = time.strftime("%Y-%m-%dT%H:%M:%S", time.gmtime(t_sub + MAX_DUR_S / 2))
            print(f"[pictor-bot]   analyze later:  hi-spectrum <csv> --pictor "
                  f"--ra {ra_mid:.2f} --dec 38.0 --time {t_mid_iso}")
    print("\n[pictor-bot] all done -- check your inbox in the morning")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
