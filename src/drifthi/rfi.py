"""RFI flagging on the calibrated (time, channel) quotient spectra.

RFI is narrowband and fixed in topocentric frequency, so flagging happens in
channel space before any Doppler regridding. Three passes:
  1. whole cycles with wild broadband power (lawnmowers, arcing, gain jumps)
  2. channels whose time-variance is far above the radiometer expectation
     (persistent carriers)
  3. individual (cycle, channel) outliers (bursty interference)

Channels inside `protect` (the velocity window where galactic HI can live)
are exempt from passes 2 and 3 -- the line is a genuine many-sigma outlier
and would otherwise flag itself.
"""

from __future__ import annotations

import numpy as np


def _mad_std(x, axis=None):
    med = np.nanmedian(x, axis=axis, keepdims=True)
    return 1.4826 * np.nanmedian(np.abs(x - med), axis=axis), med


def flag_rfi(q: np.ndarray, good_chan: np.ndarray, protect: np.ndarray,
             cycle_sigma: float = 5.0, chan_sigma: float = 6.0,
             pixel_sigma: float = 6.0):
    """Return (bad_mask (N, nfft), bad_cycles (N,))."""
    n, nfft = q.shape
    bad = np.zeros((n, nfft), dtype=bool)
    bad[:, ~good_chan] = True
    ref = good_chan & ~protect          # line-free reference channels

    # 1. cycles with deviant broadband power (judged on line-free channels)
    tot = np.nanmedian(np.where(ref[None, :], q, np.nan), axis=1)
    s, med = _mad_std(tot)
    bad_cycles = np.abs(tot - med.squeeze()) > cycle_sigma * max(float(s), 1e-12)

    # 2. persistently noisy channels (time-std per channel vs typical)
    ok = ~bad_cycles
    sample = q[ok] if ok.sum() >= 8 else q
    chan_std, _ = _mad_std(sample, axis=0)
    s2, med2 = _mad_std(chan_std[ref])
    noisy = chan_std > (med2 + chan_sigma * max(float(s2), 1e-12))
    bad[:, noisy & ~protect] = True

    # 3. per-pixel outliers against a per-channel median level
    med_c = np.nanmedian(sample, axis=0)
    resid = q - med_c[None, :]
    s3, _ = _mad_std(resid[:, ref])
    hits = np.abs(resid) > pixel_sigma * max(float(s3), 1e-12)
    hits[:, protect] = False
    bad |= hits

    frac = bad[:, good_chan].mean() if good_chan.any() else 0.0
    print(f"[rfi] flagged {100*frac:.2f}% of usable pixels, "
          f"{int(bad_cycles.sum())}/{n} cycles dropped")
    return bad, bad_cycles
