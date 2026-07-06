"""Hardware bench check -- verifies the SDR + LNA chain with no dish needed.

Runs against a live rtl_tcp and reports:
  1. ADC health: sample mean/spread, clipping at the rails
  2. bias-tee test: broadband power with LNA power off vs on -- a SAWbird
     that wakes up adds several dB of band-limited noise, so this proves the
     LNA is powered and alive
  3. bandpass spectra at the ON and OFF tunings + the calibration quotient,
     with the HI band marked, saved as PNG

Usage:
    tools\\rtl-sdr\\x64\\rtl_tcp.exe -a 127.0.0.1 -p 1234    (Windows)
    rtl_tcp -a 127.0.0.1 -p 1234                             (Pi)
    hi-check --config config.yaml
"""

from __future__ import annotations

import argparse
import json
import pathlib
import time

import numpy as np

from .config import load_config
from .sdr_tcp import connect_autostart


def _read_raw(sdr: RtlTcp, n_samples: int) -> np.ndarray:
    return np.frombuffer(sdr.read_samples_raw(n_samples * 2), dtype=np.uint8)


def _spectrum(sdr: RtlTcp, seconds: float, fs: float, nfft: int) -> np.ndarray:
    window = np.hanning(nfft).astype(np.float32)
    n_spec = max(8, int(seconds * fs) // nfft)
    acc = np.zeros(nfft)
    done = 0
    per_read = max(1, (1 << 20) // (2 * nfft))
    while done < n_spec:
        m = min(per_read, n_spec - done)
        d = _read_raw(sdr, m * nfft).astype(np.float32) - 127.5
        iq = (d[0::2] + 1j * d[1::2]).astype(np.complex64)
        s = np.fft.fft(iq.reshape(m, nfft) * window, axis=1)
        acc += np.sum(s.real**2 + s.imag**2, axis=0)
        done += m
    return np.fft.fftshift(acc / done)


def _broadband_power(sdr: RtlTcp, seconds: float, fs: float) -> float:
    d = _read_raw(sdr, int(seconds * fs)).astype(np.float32) - 127.5
    return float(np.mean(d**2))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="SDR + LNA bench check (no dish needed)")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--host", default=None)
    ap.add_argument("--port", type=int, default=None)
    ap.add_argument("--out", default="checks")
    ap.add_argument("--bias-test", action="store_true",
                    help="force the bias-tee on/off test even if sdr.bias_tee "
                         "is false in the config")
    ap.add_argument("--skip-bias-test", action="store_true",
                    help="don't toggle the bias tee (e.g. LNA powered externally)")
    args = ap.parse_args(argv)

    cfg = load_config(args.config)
    s = cfg["sdr"]
    host = args.host or s["host"]
    port = args.port or s["port"]
    fs, nfft = float(s["sample_rate_hz"]), int(s["nfft"])
    f_on, f_off = float(s["freq_on_hz"]), float(s["freq_off_hz"])

    out = pathlib.Path(args.out) / time.strftime("%Y%m%d_%H%M%S", time.gmtime())
    out.mkdir(parents=True, exist_ok=True)
    report: dict = {"t_unix": time.time(), "gain_db": float(s["gain_db"])}
    verdicts: list[str] = []

    print(f"[check] connecting to rtl_tcp at {host}:{port} ...")
    sdr, rtl_proc = connect_autostart(host, port)
    print(f"[check] connected: magic={sdr.magic!r}, tuner type={sdr.tuner_type}, "
          f"{sdr.tuner_gain_count} gain steps")
    report["tuner_type"] = sdr.tuner_type
    if sdr.magic != "RTL0":
        verdicts.append("FAIL: bad protocol magic -- is that really rtl_tcp?")

    sdr.set_sample_rate(fs)
    sdr.set_gain_db(float(s["gain_db"]))
    sdr.set_ppm(int(s.get("ppm", 0)))
    # bring the chain into its operating state BEFORE judging levels
    if bool(s.get("bias_tee", False)):
        sdr.set_bias_tee(True)
        time.sleep(0.5)
    sdr.set_freq(f_on)
    sdr.flush(int(0.3 * fs) * 2)

    # ---- 1. bias-tee / LNA aliveness test --------------------------------------
    # never inject DC up the coax unless the config says the LNA is bias-tee
    # powered (or the user forces it): externally powered LNAs keep it off
    do_bias = args.bias_test or (bool(s.get("bias_tee", False))
                                 and not args.skip_bias_test)
    if do_bias:
        print("[check] bias-tee test: measuring power with LNA off, then on ...")
        sdr.set_bias_tee(False)
        time.sleep(0.5)
        sdr.flush(int(0.3 * fs) * 2)
        p0 = _broadband_power(sdr, 0.5, fs)
        sdr.set_bias_tee(True)
        time.sleep(0.5)
        sdr.flush(int(0.3 * fs) * 2)
        p1 = _broadband_power(sdr, 0.5, fs)
        jump_db = 10 * np.log10(max(p1, 1e-12) / max(p0, 1e-12))
        report["bias_tee_jump_db"] = float(jump_db)
        print(f"[check] power bias-off -> bias-on: {jump_db:+.2f} dB")
        if jump_db > 1.0:
            verdicts.append(f"PASS: LNA responds to bias tee ({jump_db:+.1f} dB)")
        elif jump_db > 0.3:
            verdicts.append(f"WARN: small bias-tee power jump ({jump_db:+.2f} dB) -- "
                            "LNA may be powered but check connections")
        else:
            verdicts.append(f"FAIL?: no power change with bias tee ({jump_db:+.2f} dB) -- "
                            "SAWbird not powered, not connected, or powered externally")
    else:
        print("[check] bias-tee test skipped (sdr.bias_tee is false -- "
              "LNA powered externally)")
        sdr.set_bias_tee(False)
        verdicts.append("INFO: bias tee kept OFF; verify the LNA's power LED "
                        "is lit from its external supply")

    # ---- 2. ADC health (in the final operating state, LNA powered) -------------
    sdr.flush(int(0.3 * fs) * 2)
    raw = _read_raw(sdr, int(0.5 * fs))
    mean, std = float(raw.mean()), float(raw.std())
    rails = float(np.mean((raw <= 1) | (raw >= 254)) * 100)
    report["adc"] = {"mean": mean, "std": std, "rail_pct": rails}
    print(f"[check] ADC: mean={mean:.1f} (ideal ~127.5), std={std:.1f}, "
          f"samples at rails: {rails:.3f}%")
    if rails > 1.0:
        verdicts.append(f"FAIL: ADC clipping ({rails:.1f}% at rails) -- lower sdr.gain_db")
    elif rails > 0.05:
        verdicts.append(f"WARN: some ADC clipping ({rails:.2f}%) -- consider lower gain")
    if std < 3.0:
        verdicts.append(f"WARN: very low ADC drive (std={std:.1f}) -- raise sdr.gain_db "
                        "or check the LNA is powered")
    if not (110 < mean < 145):
        verdicts.append(f"WARN: ADC mean {mean:.1f} far from 127.5 (DC offset?)")

    # ---- 3. ON/OFF spectra + quotient -------------------------------------------
    print("[check] taking 4 s spectra at ON and OFF tunings ...")
    sdr.set_freq(f_on)
    sdr.flush(int(0.3 * fs) * 2)
    sp_on = _spectrum(sdr, 4.0, fs, nfft)
    sdr.set_freq(f_off)
    sdr.flush(int(0.3 * fs) * 2)
    sp_off = _spectrum(sdr, 4.0, fs, nfft)
    sdr.close()
    if rtl_proc is not None:
        rtl_proc.terminate()
        print("[check] stopped the rtl_tcp we started")

    from scipy.signal import savgol_filter
    off_sm = savgol_filter(sp_off, min(129, (nfft // 16) | 1), 3)
    off_sm = np.maximum(off_sm, 1e-3 * np.median(sp_off))
    q = (sp_on - sp_off) / off_sm
    inner = slice(int(0.1 * nfft), int(0.9 * nfft))
    qi = q[inner]
    q_rms = float(1.4826 * np.median(np.abs(qi - np.median(qi))))
    ripple_db = float(10 * np.log10(np.percentile(sp_on[inner], 95)
                                    / np.percentile(sp_on[inner], 5)))
    # radiometer expectation for this integration: ON and OFF each contribute
    # 1/sqrt(B*t); Hann window widens the noise bandwidth by ~1.5x
    b_chan = fs / nfft * 1.5
    q_expect = float(np.sqrt(2.0 / (b_chan * 4.0)))
    report["quotient_rms"] = q_rms
    report["quotient_rms_expected"] = q_expect
    report["bandpass_ripple_db"] = ripple_db
    print(f"[check] bandpass ripple (5-95%): {ripple_db:.1f} dB, "
          f"quotient rms: {q_rms:.2e} (radiometer expectation {q_expect:.2e})")
    ratio = q_rms / q_expect
    if 0.7 < ratio < 2.0:
        verdicts.append(f"PASS: quotient noise is at the radiometer limit "
                        f"(x{ratio:.2f} of ideal) -- the whole chain works")
    elif ratio >= 2.0:
        verdicts.append(f"WARN: quotient noise {ratio:.1f}x above the radiometer "
                        "limit -- RFI or gain instability")
    else:
        verdicts.append(f"WARN: quotient noise suspiciously low (x{ratio:.2f})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    freqs = (f_on + np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fs))) / 1e6
    freqs_off = (f_off + np.fft.fftshift(np.fft.fftfreq(nfft, 1 / fs))) / 1e6
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(10, 8))
    a1.plot(freqs, 10 * np.log10(np.maximum(sp_on, 1e-20)), lw=0.7,
            label=f"ON  @ {f_on/1e6:.1f} MHz")
    a1.plot(freqs_off, 10 * np.log10(np.maximum(sp_off, 1e-20)), lw=0.7,
            label=f"OFF @ {f_off/1e6:.1f} MHz")
    a1.axvline(1420.40575, color="r", ls="--", lw=0.8, label="HI rest")
    a1.set_xlabel("frequency [MHz]")
    a1.set_ylabel("power [dB, arb]")
    a1.set_title("raw bandpasses (should be smooth humps; carriers = RFI)")
    a1.legend(fontsize=8)
    a1.grid(alpha=0.3)
    a2.plot(freqs, q, lw=0.7, color="k")
    a2.axvspan(1420.40575 - 0.75, 1420.40575 + 0.75, color="r", alpha=0.12,
               label="galactic HI can appear here")
    a2.set_xlabel("frequency [MHz]")
    a2.set_ylabel("(ON-OFF)/OFF")
    a2.set_title("calibration quotient (flat + noise without a dish)")
    a2.legend(fontsize=8)
    a2.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "bench_spectra.png", dpi=130)
    plt.close(fig)

    report["verdicts"] = verdicts
    with open(out / "report.json", "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2)

    print("\n[check] ================= VERDICTS =================")
    for v in verdicts:
        print(f"[check] {v}")
    print(f"[check] plots + report in {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
