"""Continuous drift-scan observer. Runs on the Raspberry Pi.

Talks to a local rtl_tcp instance, integrates FFT power spectra, and
frequency-switches between an ON band (containing the 1420.406 MHz HI line)
and a line-free OFF band. The OFF spectrum is used later to divide out the
receiver bandpass, which is what makes HI visible with an RTL-SDR at all.

Start the driver first (rtl-sdr-blog fork for a V4 dongle):
    rtl_tcp -a 127.0.0.1 -p 1234
then:
    hi-observe --config config.yaml

Only needs numpy + pyyaml on the Pi. Ctrl+C stops cleanly; driver hiccups
trigger an automatic reconnect, so it can run unattended for days.
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import sys
import time

import numpy as np

from . import __version__
from .config import load_config
from .sdr_tcp import RtlTcp, connect_autostart
from . import session as sess


def _integrate(sdr: RtlTcp, n_samples: int, nfft: int, window: np.ndarray):
    """Read n_samples I/Q samples and return (summed |FFT|^2, n_ffts)."""
    n_spec_total = n_samples // nfft
    acc = np.zeros(nfft, dtype=np.float64)
    n_done = 0
    # ~0.35 s of data per read at 2.4 MS/s
    spec_per_read = max(1, (1 << 20) // (2 * nfft))
    while n_done < n_spec_total:
        m = min(spec_per_read, n_spec_total - n_done)
        raw = sdr.read_samples_raw(m * nfft * 2)
        d = np.frombuffer(raw, dtype=np.uint8).astype(np.float32)
        d -= 127.5
        iq = (d[0::2] + 1j * d[1::2]).astype(np.complex64)
        blocks = iq.reshape(m, nfft) * window
        spec = np.fft.fft(blocks, axis=1)
        acc += np.sum(spec.real**2 + spec.imag**2, axis=0)
        n_done += m
    return np.fft.fftshift(acc / max(n_done, 1)), n_done


def _connect(cfg_sdr: dict, bias_tee: bool):
    sdr, proc = connect_autostart(cfg_sdr["host"], cfg_sdr["port"])
    print(f"[observe] connected to rtl_tcp at {sdr.host}:{sdr.port} "
          f"(tuner type {sdr.tuner_type}, {sdr.tuner_gain_count} gain steps)")
    sdr.set_sample_rate(float(cfg_sdr["sample_rate_hz"]))
    sdr.set_gain_db(float(cfg_sdr["gain_db"]))
    sdr.set_ppm(int(cfg_sdr.get("ppm", 0)))
    if bias_tee:
        sdr.set_bias_tee(True)
        print("[observe] bias tee ON (powering the LNA through the coax)")
    else:
        sdr.set_bias_tee(False)
        print("[observe] bias tee OFF (LNA powered externally)")
    return sdr, proc


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Continuous frequency-switched HI drift scan")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--out", default=None,
                    help="root directory for sessions (default: paths.raw_dir "
                         "from the config, else data/raw)")
    ap.add_argument("--tag", default="obs", help="session name suffix")
    ap.add_argument("--duration-h", type=float, default=0.0,
                    help="stop after this many hours (0 = run until Ctrl+C)")
    ap.add_argument("--no-bias-tee", action="store_true")
    args = ap.parse_args(argv)

    # line-buffer stdout so `journalctl -f` / log files show cycles live
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except AttributeError:
        pass

    # systemd stops services with SIGTERM: convert it to KeyboardInterrupt so
    # the same clean-shutdown path runs (flush chunk, finalize meta.json)
    def _sigterm(_sig, _frame):
        raise KeyboardInterrupt
    try:
        signal.signal(signal.SIGTERM, _sigterm)
    except (ValueError, OSError):
        pass

    cfg = load_config(args.config)
    out_root = args.out or cfg.get("paths", {}).get("raw_dir", "data/raw")
    s = cfg["sdr"]
    fs = float(s["sample_rate_hz"])
    nfft = int(s["nfft"])
    f_on = float(s["freq_on_hz"])
    f_off = float(s["freq_off_hz"])
    t_on = float(s["t_on_s"])
    t_off = float(s["t_off_s"])
    settle = float(s["settle_s"])
    chunk_cycles = int(s["chunk_cycles"])
    bias_tee = bool(s.get("bias_tee", False)) and not args.no_bias_tee

    n_on = int(t_on * fs) // nfft * nfft
    n_off = int(t_off * fs) // nfft * nfft
    n_settle = max(nfft, int(settle * fs) // nfft * nfft) * 2  # bytes = samples*2
    window = np.hanning(nfft).astype(np.float32)

    # safety: if raw_dir lives on an external mount (/mnt/.. or /media/..)
    # that isn't actually mounted, data would silently land on the SD card.
    ap_root = os.path.abspath(out_root)
    for base in ("/mnt", "/media"):
        if ap_root.startswith(base + os.sep):
            mountdir = os.sep.join(ap_root.split(os.sep)[:3])  # e.g. /mnt/ssd
            if not os.path.ismount(mountdir):
                print(f"[observe] WARNING: {mountdir} is NOT a mounted filesystem "
                      f"-- data would go to the SD card, not the drive. "
                      f"Run 'sudo mount -a' (or plug the drive in before boot), "
                      f"then restart.", file=sys.stderr)
            break

    out = sess.new_session_dir(out_root, args.tag)
    meta = {
        "version": __version__,
        "kind": "observation",
        "t_start_unix": time.time(),
        "sample_rate_hz": fs,
        "freq_on_hz": f_on,
        "freq_off_hz": f_off,
        "nfft": nfft,
        "gain_db": float(s["gain_db"]),
        "t_on_s": t_on,
        "t_off_s": t_off,
        "site": cfg["site"],
        "pointing": cfg["pointing"],
        "config": cfg,
    }
    sess.write_meta(out, meta)
    print(f"[observe] session dir: {out}")
    print(f"[observe] ON {f_on/1e6:.3f} MHz / OFF {f_off/1e6:.3f} MHz, "
          f"{t_on:.0f}s + {t_off:.0f}s cycles, nfft={nfft} "
          f"({fs/nfft:.0f} Hz = {299792.458*fs/nfft/1.42040575e9:.3f} km/s per channel)")

    deadline = time.time() + args.duration_h * 3600 if args.duration_h > 0 else None
    cyc_t, cyc_on, cyc_off, cyc_non, cyc_noff = [], [], [], [], []
    chunk_idx = 0
    cycle = 0
    sdr = None
    rtl_proc = None   # rtl_tcp we spawned ourselves (terminated on exit)
    p_ref = None      # first cycle's power: reference for the live dB column

    def flush():
        nonlocal chunk_idx
        if cyc_t:
            sess.write_chunk(out, chunk_idx, cyc_t, cyc_on, cyc_off, cyc_non, cyc_noff)
            print(f"[observe] wrote chunk {chunk_idx:06d} ({len(cyc_t)} cycles)")
            chunk_idx += 1
            for lst in (cyc_t, cyc_on, cyc_off, cyc_non, cyc_noff):
                lst.clear()

    try:
        while deadline is None or time.time() < deadline:
            try:
                if sdr is None:
                    sdr, proc = _connect(s, bias_tee)
                    rtl_proc = proc or rtl_proc
                t0 = time.time()
                sdr.set_freq(f_on)
                sdr.flush(n_settle)
                spec_on, m_on = _integrate(sdr, n_on, nfft, window)
                sdr.set_freq(f_off)
                sdr.flush(n_settle)
                spec_off, m_off = _integrate(sdr, n_off, nfft, window)
                t1 = time.time()

                cyc_t.append(0.5 * (t0 + t1))
                cyc_on.append(spec_on.astype(np.float32))
                cyc_off.append(spec_off.astype(np.float32))
                cyc_non.append(m_on)
                cyc_noff.append(m_off)
                cycle += 1

                ratio = float(np.median(spec_on) / np.median(spec_off))
                p_now = float(np.median(spec_off))
                if p_ref is None:
                    p_ref = p_now
                p_db = 10 * np.log10(max(p_now, 1e-30) / p_ref)
                print(f"[observe] cycle {cycle:6d}  {time.strftime('%H:%M:%S', time.gmtime(t1))}Z  "
                      f"P_on/P_off={ratio:.4f}  P={p_db:+5.2f} dB  ({t1-t0:.1f}s)")
                if len(cyc_t) >= chunk_cycles:
                    flush()
            except (ConnectionError, OSError, TimeoutError) as e:
                print(f"[observe] SDR error: {e!r} -- reconnecting in 5 s", file=sys.stderr)
                if sdr is not None:
                    sdr.close()
                    sdr = None
                flush()
                time.sleep(5)
    except KeyboardInterrupt:
        print("\n[observe] stopped by user")
    finally:
        flush()
        if sdr is not None:
            sdr.close()
        if rtl_proc is not None:
            rtl_proc.terminate()
            print("[observe] stopped the rtl_tcp we started")
        meta["t_end_unix"] = time.time()
        meta["n_cycles"] = cycle
        sess.write_meta(out, meta)
        print(f"[observe] done: {cycle} cycles in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
