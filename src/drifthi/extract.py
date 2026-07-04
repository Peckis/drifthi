"""Science extraction from the RA-binned spectra.

Per RA bin:
  * noise level, peak T_B and its velocity
  * emission mask, moment 0 (integrated intensity), intensity-weighted
    velocity (moment 1) and dispersion (moment 2)
  * HI column density, optically thin: N_HI = 1.823e18 * integral(T_B dv)
  * multi-Gaussian decomposition (up to 4 components)
Globally:
  * kinematic distances (flat rotation curve) -> face-on Milky Way HI map
  * tangent-point rotation curve where the sightline passes the inner Galaxy
"""

from __future__ import annotations

import csv
import pathlib

import numpy as np
from scipy.optimize import curve_fit

R0_KPC = 8.178   # GRAVITY 2019
V0_KMS = 220.0   # flat rotation speed


def noise_sigma(spec: np.ndarray, v: np.ndarray, v_quiet: float = 180.0) -> float:
    quiet = np.abs(v) > v_quiet
    x = spec[quiet & np.isfinite(spec)]
    if x.size < 10:
        x = spec[np.isfinite(spec)]
    if x.size < 10:
        return np.nan
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def emission_mask(spec: np.ndarray, sigma: float, thresh: float = 3.0,
                  grow: int = 3) -> np.ndarray:
    m = np.nan_to_num(spec) > thresh * sigma
    if grow > 0 and m.any():
        idx = np.flatnonzero(m)
        for i in idx:
            m[max(0, i - grow): i + grow + 1] = True
    return m


def _multi_gauss(v, *p):
    y = np.zeros_like(v)
    for i in range(0, len(p), 3):
        a, v0, s = p[i], p[i + 1], max(p[i + 2], 0.3)
        y = y + a * np.exp(-0.5 * ((v - v0) / s) ** 2)
    return y


def decompose(v: np.ndarray, spec: np.ndarray, sigma: float,
              max_comp: int = 4, min_snr: float = 5.0) -> list[dict]:
    """Iterative peak-find + joint refit. Returns list of {amp, v0, sigma_v}."""
    from scipy.ndimage import uniform_filter1d

    ok = np.isfinite(spec)
    if ok.sum() < 20 or not np.isfinite(sigma) or sigma <= 0:
        return []
    vv, ss = v[ok], spec[ok].copy()
    resid = ss.copy()
    p0 = []
    for _ in range(max_comp):
        smooth = uniform_filter1d(resid, 5)
        i = int(np.argmax(smooth))
        if smooth[i] < min_snr * sigma:
            break
        # initial width from half-max extent around the peak
        half = smooth[i] / 2
        l = i
        while l > 0 and smooth[l] > half:
            l -= 1
        r = i
        while r < len(vv) - 1 and smooth[r] > half:
            r += 1
        w = max((vv[r] - vv[l]) / 2.3548, 1.5)
        p0 += [smooth[i], vv[i], w]
        resid = ss - _multi_gauss(vv, *p0)
    if not p0:
        return []
    try:
        n = len(p0) // 3
        lo = [0.0, vv.min(), 0.5] * n
        hi = [np.inf, vv.max(), 80.0] * n
        popt, _ = curve_fit(_multi_gauss, vv, ss, p0=p0, bounds=(lo, hi), maxfev=8000)
    except (RuntimeError, ValueError):
        popt = np.array(p0)
    comps = [{"amp_k": float(popt[i]), "v0_kms": float(popt[i + 1]),
              "sigma_v_kms": float(abs(popt[i + 2]))}
             for i in range(0, len(popt), 3)]
    return sorted(comps, key=lambda c: -c["amp_k"])


def analyze_bins(v: np.ndarray, W: np.ndarray, ra_bins: np.ndarray,
                 dec_bins: np.ndarray, l_bins: np.ndarray, b_bins: np.ndarray,
                 out_csv: pathlib.Path) -> list[dict]:
    dv = float(np.median(np.diff(v)))
    rows = []
    for i in range(W.shape[0]):
        spec = W[i]
        if not np.isfinite(spec).any():
            continue
        sig = noise_sigma(spec, v)
        m = emission_mask(spec, sig) if np.isfinite(sig) else np.zeros_like(spec, bool)
        r = {"ra_deg": float(ra_bins[i]), "dec_deg": float(dec_bins[i]),
             "l_deg": float(l_bins[i]), "b_deg": float(b_bins[i]),
             "noise_k": sig, "peak_k": float(np.nanmax(spec)),
             "v_peak_kms": float(v[np.nanargmax(spec)])}
        if m.any():
            sm = np.nan_to_num(spec)
            mom0 = float(np.sum(sm[m]) * dv)                       # K km/s
            mom1 = float(np.sum(sm[m] * v[m]) / max(np.sum(sm[m]), 1e-12))
            mom2 = float(np.sqrt(max(np.sum(sm[m] * (v[m] - mom1) ** 2)
                                     / max(np.sum(sm[m]), 1e-12), 0.0)))
            r.update(mom0_k_kms=mom0, v_centroid_kms=mom1, v_disp_kms=mom2,
                     nhi_cm2=1.823e18 * mom0,
                     v_min_emis_kms=float(v[m].min()), v_max_emis_kms=float(v[m].max()))
        else:
            r.update(mom0_k_kms=0.0, v_centroid_kms=np.nan, v_disp_kms=np.nan,
                     nhi_cm2=0.0, v_min_emis_kms=np.nan, v_max_emis_kms=np.nan)
        comps = decompose(v, spec, sig)
        for j in range(4):
            if j < len(comps):
                r[f"g{j+1}_amp_k"] = comps[j]["amp_k"]
                r[f"g{j+1}_v0_kms"] = comps[j]["v0_kms"]
                r[f"g{j+1}_sigma_kms"] = comps[j]["sigma_v_kms"]
            else:
                r[f"g{j+1}_amp_k"] = r[f"g{j+1}_v0_kms"] = r[f"g{j+1}_sigma_kms"] = np.nan
        rows.append(r)

    if rows:
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            wcsv = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            wcsv.writeheader()
            wcsv.writerows(rows)
        print(f"[extract] wrote {out_csv} ({len(rows)} RA bins)")
    return rows


def faceon_points(v: np.ndarray, W: np.ndarray, l_bins: np.ndarray,
                  b_bins: np.ndarray, snr_thresh: float = 4.0):
    """Kinematic (x, y, T_B) points in the galactic plane, outer Galaxy only.

    Flat rotation: v_lsr = V0 sin(l) (R0/R - 1)  =>  R = R0 / (1 + v/(V0 sin l)).
    Outer-Galaxy solutions (R > R0) have a unique distance:
        d = R0 cos(l) + sqrt(R^2 - R0^2 sin^2 l).
    Sun at (0, R0); Galactic center at origin; x = d sin l, y = R0 - d cos l.
    """
    xs, ys, ts = [], [], []
    for i in range(W.shape[0]):
        spec = W[i]
        sig = noise_sigma(spec, v)
        if not np.isfinite(sig) or abs(b_bins[i]) > 20:
            continue
        l = np.radians(l_bins[i])
        sinl = np.sin(l)
        if abs(sinl) < 0.15:
            continue
        good = np.nan_to_num(spec) > snr_thresh * sig
        for k in np.flatnonzero(good):
            ratio = 1.0 + v[k] / (V0_KMS * sinl)
            if ratio <= 0:
                continue
            R = R0_KPC / ratio
            if R <= R0_KPC * 1.02 or R > 30:
                continue  # keep only unambiguous outer-Galaxy gas
            disc = R * R - (R0_KPC * sinl) ** 2
            if disc < 0:
                continue
            d = R0_KPC * np.cos(l) + np.sqrt(disc)
            if d <= 0 or d > 35:
                continue
            xs.append(d * np.sin(l))
            ys.append(R0_KPC - d * np.cos(l))
            ts.append(spec[k])
    return np.array(xs), np.array(ys), np.array(ts)


def rotation_curve(v: np.ndarray, W: np.ndarray, l_bins: np.ndarray,
                   b_bins: np.ndarray):
    """Tangent-point V(R) for sightlines through the inner Galaxy (|b|<10)."""
    out = []
    for i in range(W.shape[0]):
        l_deg = l_bins[i] % 360.0
        if abs(b_bins[i]) > 10:
            continue
        spec = W[i]
        sig = noise_sigma(spec, v)
        if not np.isfinite(sig):
            continue
        good = np.nan_to_num(spec) > 4.0 * sig
        if not good.any():
            continue
        l = np.radians(l_deg)
        if 0 < l_deg < 90:
            v_t = float(v[good].max())
            vc = v_t + V0_KMS * np.sin(l)
        elif 270 < l_deg < 360:
            v_t = float(v[good].min())
            vc = -(v_t + V0_KMS * np.sin(l))  # sin(l)<0; terminal is most negative
        else:
            continue
        R = R0_KPC * abs(np.sin(l))
        if 0.5 < R < R0_KPC and 30 < vc < 350:
            out.append((R, vc, l_deg))
    return np.array(out) if out else np.empty((0, 3))
