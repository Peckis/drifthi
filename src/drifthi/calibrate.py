"""Frequency-switched calibration and baseline removal.

The quotient (ON - OFF) / smooth(OFF) cancels, per channel:
  * the receiver bandpass shape (multiplicative, identical relative to the LO
    at both tunings), and
  * additive artifacts present at both tunings (DC spike, ADC spurs).
The result is T_line / T_sys, converted to Kelvin with an assumed T_sys
(later refined by hi-compare against HI4PI).
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import median_filter
from scipy.signal import savgol_filter


def channel_mask(nfft: int, edge_frac: float, dc_halfwidth: int) -> np.ndarray:
    """True where a channel is usable (inside filter rolloff, away from DC)."""
    good = np.ones(nfft, dtype=bool)
    ne = int(edge_frac * nfft)
    if ne > 0:
        good[:ne] = False
        good[-ne:] = False
    c = nfft // 2
    good[c - dc_halfwidth: c + dc_halfwidth + 1] = False
    return good


def quotient_kelvin(on: np.ndarray, off: np.ndarray, tsys_k: float) -> np.ndarray:
    """(N, nfft) ON/OFF power -> antenna temperature estimate in K."""
    nfft = on.shape[1]
    win = min(129, (nfft // 16) | 1)
    off_smooth = savgol_filter(median_filter(off, size=(1, 9)), win, 3, axis=1)
    # savgol can overshoot <=0 on the steep rolloff shoulders; floor at a
    # fraction of each spectrum's median so the quotient can't explode there
    floor = 1e-3 * np.median(off, axis=1, keepdims=True)
    off_smooth = np.maximum(off_smooth, np.maximum(floor, np.finfo(np.float32).tiny))
    return (on - off) / off_smooth * tsys_k


def fit_baseline(spec: np.ndarray, x: np.ndarray, order: int = 3,
                 n_iter: int = 5, clip_hi: float = 2.0, clip_lo: float = 3.0):
    """Sigma-clipped polynomial baseline of a single spectrum.

    Asymmetric clipping: emission is positive, so points far above the fit are
    excluded aggressively while the fit still rides on the noise floor.
    Returns (baseline, keep_mask); NaNs in `spec` are ignored.
    """
    keep = np.isfinite(spec)
    if keep.sum() < (order + 1) * 3:
        return np.zeros_like(spec), keep
    xs = (x - x.mean()) / max(float(np.ptp(x)), 1e-9)
    base = np.zeros_like(spec)
    for _ in range(n_iter):
        c = np.polynomial.polynomial.polyfit(xs[keep], spec[keep], order)
        base = np.polynomial.polynomial.polyval(xs, c)
        r = spec - base
        s = 1.4826 * np.nanmedian(np.abs(r[keep] - np.nanmedian(r[keep])))
        if s <= 0 or not np.isfinite(s):
            break
        new = np.isfinite(spec) & (r < clip_hi * s) & (r > -clip_lo * s)
        if new.sum() < (order + 1) * 3 or np.array_equal(new, keep):
            keep = new if new.sum() >= (order + 1) * 3 else keep
            break
        keep = new
    return base, keep


def remove_baselines(T: np.ndarray, v: np.ndarray, order: int = 3) -> np.ndarray:
    out = np.empty_like(T)
    for i in range(T.shape[0]):
        base, _ = fit_baseline(T[i], v, order=order)
        out[i] = T[i] - base
    return out
