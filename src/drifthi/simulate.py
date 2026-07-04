"""Generate a synthetic raw session with realistic instrument effects.

Lets you test hi-process and hi-compare end-to-end before the telescope is
finished. The sky comes from the HI4PI strip cache if present, otherwise a
built-in toy Milky Way (galactic plane + local gas + a few clouds), so it
works fully offline.

Instrument model: receiver bandpass with ripple and edge rolloff, slow gain
drift, DC spike, radiometer noise for the actual integration time, a
persistent RFI carrier and random RFI bursts. Optionally inject a "true"
pointing error (--true-daz/--true-del) that hi-compare should recover.

Usage:
    hi-simulate --duration-h 24 --true-daz 5 --true-del -4
"""

from __future__ import annotations

import argparse
import pathlib
import time

import numpy as np

from . import __version__, hi4pi, velocity
from . import session as sess
from .config import load_config


def toy_strip(cfg: dict, dec_min: float, dec_max: float,
              out_path: pathlib.Path) -> pathlib.Path:
    """Analytic HI sky in the same strip-cache format as hi4pi.build_strip_cache."""
    h = cfg["hi4pi"]
    sky_bin = float(h["sky_bin_deg"])
    dv = 1.28821497 * int(h["v_bin_ch"])
    v = np.arange(-float(h["v_max_kms"]), float(h["v_max_kms"]) + dv / 2, dv)
    ra = np.arange(0.0, 360.0, sky_bin) + sky_bin / 2
    dec = np.arange(dec_min, dec_max, sky_bin) + sky_bin / 2

    rr, dd = np.meshgrid(ra, dec, indexing="ij")
    l, b = velocity.galactic_lb(rr.ravel(), dd.ravel())
    l = l.reshape(rr.shape)
    b = b.reshape(rr.shape)
    lr = np.radians(l)

    def gauss_v(amp, v0, sig):
        return amp[..., None] * np.exp(-0.5 * ((v[None, None, :]
                                                - v0[..., None]) / sig) ** 2)

    # galactic plane, outer-galaxy-flavored velocity field
    amp1 = 70.0 * np.exp(-0.5 * (b / 4.5) ** 2)
    v1 = -60.0 * np.sin(lr) * np.cos(np.radians(b))
    cube = gauss_v(amp1, v1, 14.0)
    # local, wide-latitude gas near v=0
    amp2 = 16.0 * np.exp(-0.5 * (b / 13.0) ** 2)
    v2 = -6.0 * np.sin(2 * lr)
    cube += gauss_v(amp2, v2, 7.0)
    # discrete clouds / HVC-ish features
    for l0, b0, rad, amp, v0, sv in [(140, 18, 7, 24, -85, 6),
                                     (200, 8, 6, 18, -35, 5),
                                     (95, -12, 8, 14, 25, 7),
                                     (110, 30, 10, 7, -130, 9)]:
        dl = (l - l0 + 180) % 360 - 180
        a = amp * np.exp(-0.5 * ((dl / rad) ** 2 + ((b - b0) / rad) ** 2))
        cube += gauss_v(a, np.full_like(a, float(v0)), sv)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out_path, ra=ra.astype(np.float32),
                        dec=dec.astype(np.float32), v=v.astype(np.float32),
                        cube=cube.astype(np.float32))
    print(f"[simulate] toy sky strip written: {out_path}")
    return out_path


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Simulate a raw drift-scan session")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--out", default=None,
                    help="output root (default: paths.raw_dir from config)")
    ap.add_argument("--tag", default="sim")
    ap.add_argument("--duration-h", type=float, default=24.0)
    ap.add_argument("--start-unix", type=float, default=None)
    ap.add_argument("--true-daz", type=float, default=0.0,
                    help="injected azimuth pointing error [deg]")
    ap.add_argument("--true-del", type=float, default=0.0,
                    help="injected elevation pointing error [deg]")
    ap.add_argument("--strip", default=None, help="sky strip cache to observe")
    ap.add_argument("--rfi-prob", type=float, default=0.02)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    rng = np.random.default_rng(args.seed)
    s = cfg["sdr"]
    fs, nfft = float(s["sample_rate_hz"]), int(s["nfft"])
    t_on, t_off = float(s["t_on_s"]), float(s["t_off_s"])
    tsys = float(cfg["processing"]["tsys_assumed_k"])

    # sky strip: user-supplied > HI4PI cache > toy model
    cache = pathlib.Path(cfg["hi4pi"]["cache_dir"]) / "strip_cache.npz"
    if args.strip:
        strip_path = pathlib.Path(args.strip)
    elif cache.exists():
        strip_path = cache
        print(f"[simulate] using HI4PI strip cache {cache}")
    else:
        t_probe = args.start_unix or time.time()
        _, dec0 = velocity.pointing_radec(np.array([t_probe]), cfg,
                                          args.true_daz, args.true_del)
        margin = cfg["pointing"]["beam_fwhm_deg"] + cfg["compare"]["search_deg"] + 6
        strip_path = pathlib.Path(cfg["hi4pi"]["cache_dir"]) / "toy_strip.npz"
        toy_strip(cfg, float(dec0[0]) - margin, float(dec0[0]) + margin, strip_path)

    ra_g, dec_g, v_g, cube = hi4pi.load_strip(strip_path)
    sm = hi4pi.beam_smooth_strip(ra_g, dec_g, cube,
                                 float(cfg["pointing"]["beam_fwhm_deg"]))

    # cycle timing and true (injected) pointing track
    period = t_on + t_off + 2 * float(s["settle_s"]) + 0.6
    n_cyc = max(4, int(args.duration_h * 3600 / period))
    t0 = args.start_unix or time.time()
    t_mid = t0 + period * (np.arange(n_cyc) + 0.5)
    ra_t, dec_t = velocity.pointing_radec(t_mid, cfg, args.true_daz, args.true_del)
    vcorr = velocity.vlsr_correction_kms(t_mid, ra_t, dec_t, cfg)
    Ta_track = hi4pi.track_spectra(ra_g, dec_g, v_g, sm, ra_t, dec_t, v_g)
    Ta_track = np.nan_to_num(Ta_track, nan=0.0)
    print(f"[simulate] {n_cyc} cycles over {args.duration_h:.1f} h, "
          f"true pointing error daz={args.true_daz:+.1f} del={args.true_del:+.1f} deg, "
          f"track dec ~ {np.mean(dec_t):+.1f} deg")

    # instrument model in channel space
    freqs = velocity.channel_freqs_hz(float(s["freq_on_hz"]), fs, nfft)
    v_topo = velocity.v_topo_kms(freqs)
    xn = np.linspace(-1.0, 1.0, nfft)
    edge = np.ones(nfft)
    ramp = np.clip((np.abs(xn) - 0.8) / 0.2, 0, 1)
    edge *= 0.03 + 0.97 * 0.5 * (1 + np.cos(np.pi * ramp))
    G = 1e-3 * (1 + 0.12 * np.sin(2 * np.pi * 2.5 * xn + 0.7) + 0.05 * xn) * edge
    tsys_f = tsys * (1 + 0.03 * xn)
    dc = np.zeros(nfft)
    dc[nfft // 2] = 30 * tsys * 1e-3
    dc[nfft // 2 - 1] = dc[nfft // 2 + 1] = 8 * tsys * 1e-3

    m_on = int(t_on * fs) // nfft
    m_off = int(t_off * fs) // nfft
    # persistent carrier outside the protected velocity window
    k_persist = int(np.argmin(np.abs(v_topo - 265.0)))

    out = sess.new_session_dir(
        args.out or cfg.get("paths", {}).get("raw_dir", "data/raw"), args.tag)
    meta = {
        "version": __version__, "kind": "simulation",
        "t_start_unix": float(t_mid[0]), "sample_rate_hz": fs,
        "freq_on_hz": float(s["freq_on_hz"]), "freq_off_hz": float(s["freq_off_hz"]),
        "nfft": nfft, "gain_db": float(s["gain_db"]),
        "t_on_s": t_on, "t_off_s": t_off,
        "site": cfg["site"], "pointing": cfg["pointing"],
        "true_daz_deg": args.true_daz, "true_del_deg": args.true_del,
        "strip": str(strip_path), "seed": args.seed,
    }
    sess.write_meta(out, meta)

    chunk_cycles = int(s["chunk_cycles"])
    buf = {k: [] for k in ("t", "on", "off", "non", "noff")}
    chunk_idx = 0
    for i in range(n_cyc):
        Ta_chan = np.interp(v_topo + vcorr[i], v_g, Ta_track[i], left=0.0, right=0.0)
        drift = 1 + 0.03 * np.sin(2 * np.pi * (t_mid[i] - t0) / (5 * 3600))
        p_on = drift * G * (tsys_f + Ta_chan) + dc
        p_off = drift * G * tsys_f + dc
        on = p_on * (1 + rng.standard_normal(nfft) / np.sqrt(m_on))
        off = p_off * (1 + rng.standard_normal(nfft) / np.sqrt(m_off))
        on[k_persist] *= 1 + 3.0 + 0.8 * rng.standard_normal()
        if rng.random() < args.rfi_prob:
            k = rng.integers(int(0.1 * nfft), int(0.9 * nfft))
            on[k:k + rng.integers(1, 4)] *= rng.uniform(3, 40)
        if rng.random() < args.rfi_prob:
            k = rng.integers(int(0.1 * nfft), int(0.9 * nfft))
            off[k:k + rng.integers(1, 4)] *= rng.uniform(3, 40)

        buf["t"].append(t_mid[i])
        buf["on"].append(on.astype(np.float32))
        buf["off"].append(off.astype(np.float32))
        buf["non"].append(m_on)
        buf["noff"].append(m_off)
        if len(buf["t"]) >= chunk_cycles or i == n_cyc - 1:
            sess.write_chunk(out, chunk_idx, buf["t"], buf["on"], buf["off"],
                             buf["non"], buf["noff"])
            chunk_idx += 1
            for lst in buf.values():
                lst.clear()
        if (i + 1) % 500 == 0:
            print(f"[simulate] {i+1}/{n_cyc} cycles")

    meta["t_end_unix"] = float(t_mid[-1])
    meta["n_cycles"] = n_cyc
    sess.write_meta(out, meta)
    print(f"[simulate] session written: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
