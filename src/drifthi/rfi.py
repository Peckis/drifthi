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


def flag_persistent_channels(on: np.ndarray, off: np.ndarray,
                             good_chan: np.ndarray, sigma: float = 5.0,
                             max_width: int = 80, pad: int = 8) -> np.ndarray:
    """Channels holding persistent (all-night) RFI: carriers and SMPS combs.

    Works on the RAW ON and OFF spectra, median-combined over time. RFI is
    fixed in channel space so it survives the median at full contrast, while
    the HI line drifts with the LSR correction (hundreds of channels per
    session) and smears out. Features wider than `max_width` channels are
    ignored as a safety net (a broad line remnant can't be flagged).
    Returns a bool mask over channels (True = RFI).
    """
    from scipy.signal import savgol_filter

    nfft = on.shape[1]
    bad = np.zeros(nfft, dtype=bool)
    win = min(257, (nfft // 8) | 1)
    for arr in (on, off):
        med = np.median(arr, axis=0)
        sm = np.maximum(savgol_filter(med, win, 3), 1e-30)
        r = med / sm - 1.0
        s, m = _mad_std(r[good_chan])
        hits = (r - float(m.squeeze())) > sigma * max(float(s), 1e-9)
        hits &= good_chan
        idx = np.flatnonzero(hits)
        if idx.size == 0:
            continue
        for g in np.split(idx, np.flatnonzero(np.diff(idx) > 5) + 1):
            if g[-1] - g[0] + 1 <= max_width:
                bad[max(g[0] - pad, 0): g[-1] + pad + 1] = True
    n_feat = int(np.sum(np.diff(bad.astype(int)) == 1))
    print(f"[rfi] persistent-channel flag: {int(bad.sum())} channels "
          f"in ~{n_feat} features (carriers / SMPS comb)")
    return bad


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
