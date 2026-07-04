"""Frequency <-> velocity axes, drift-scan pointing, and LSR corrections."""

from __future__ import annotations

import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

# The pipeline must run on an offline Pi; arcsecond-level IERS corrections are
# irrelevant under a 17-degree beam.
iers.conf.auto_download = False
iers.conf.auto_max_age = None

C_KMS = 299792.458
F0_HI = 1420.405751768e6  # HI rest frequency [Hz]

# Solar motion w.r.t. the kinematic LSR: 20 km/s toward 18h03m50s +30d00m (J2000)
_LSR_APEX = SkyCoord("18h03m50.29s", "+30d00m16.8s", frame="icrs")
_LSR_SPEED = 20.0  # km/s


def channel_freqs_hz(fc_hz: float, fs_hz: float, nfft: int) -> np.ndarray:
    """Sky frequency of each fftshifted channel."""
    return fc_hz + np.fft.fftshift(np.fft.fftfreq(nfft, d=1.0 / fs_hz))


def v_topo_kms(freqs_hz: np.ndarray) -> np.ndarray:
    """Topocentric radial velocity (radio convention) of each channel."""
    return C_KMS * (F0_HI - freqs_hz) / F0_HI


def earth_location(cfg: dict) -> EarthLocation:
    s = cfg["site"]
    return EarthLocation(lat=s["lat_deg"] * u.deg, lon=s["lon_deg"] * u.deg,
                         height=s.get("height_m", 0.0) * u.m)


def pointing_radec(t_unix: np.ndarray, cfg: dict,
                   daz_deg: float = 0.0, del_deg: float = 0.0):
    """ICRS (ra, dec) in degrees of the fixed az/el pointing at each time."""
    loc = earth_location(cfg)
    p = cfg["pointing"]
    t = Time(np.atleast_1d(t_unix), format="unix")
    aa = AltAz(az=(p["az_deg"] + daz_deg) * u.deg,
               alt=np.clip(p["el_deg"] + del_deg, 0.5, 89.9) * u.deg,
               obstime=t, location=loc)
    icrs = SkyCoord(aa).icrs
    return icrs.ra.deg, icrs.dec.deg


def vlsr_correction_kms(t_unix: np.ndarray, ra_deg, dec_deg, cfg: dict) -> np.ndarray:
    """Correction to ADD to topocentric velocity to get v_LSR.

    v_lsr = v_topo + barycentric correction + solar motion projected on the
    line of sight (kinematic LSR definition).
    """
    loc = earth_location(cfg)
    t = Time(np.atleast_1d(t_unix), format="unix")
    sc = SkyCoord(np.atleast_1d(ra_deg) * u.deg, np.atleast_1d(dec_deg) * u.deg)
    vbary = sc.radial_velocity_correction(kind="barycentric", obstime=t,
                                          location=loc).to_value(u.km / u.s)
    apex_term = _LSR_SPEED * np.cos(sc.separation(_LSR_APEX).rad)
    return vbary + apex_term


def galactic_lb(ra_deg, dec_deg):
    g = SkyCoord(np.atleast_1d(ra_deg) * u.deg,
                 np.atleast_1d(dec_deg) * u.deg).galactic
    return g.l.deg, g.b.deg
