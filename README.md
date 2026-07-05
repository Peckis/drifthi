# drifthi — a drift-scan HI (21 cm) pipeline for a small dish + RTL-SDR

End-to-end pipeline for a DIY hydrogen-line telescope:

> 80 cm offset TV dish → cantenna feed → SAWbird H1 (LNA + filter) →
> RTL-SDR Blog V4 → Raspberry Pi (acquisition) → any PC (processing)

No alt/az mount needed — the dish stares at a fixed azimuth/elevation and the
Earth does the scanning (drift scan). Over 24 h you sweep a full circle of
constant declination through the Milky Way.

## Tools

| command          | runs on | what it does |
|------------------|---------|--------------|
| `hi-observe`     | Pi      | continuous frequency-switched acquisition via `rtl_tcp` |
| `hi-process`     | PC      | calibration → RFI flagging → LSR regridding → baselines → waterfalls → science extraction |
| `hi-fetch-hi4pi` | PC      | downloads the HI4PI survey tiles for your declination strip, builds a compact cache |
| `hi-compare`     | PC      | cross-correlates your data with HI4PI: **fits your true pointing** and absolute calibration, makes side-by-side/overlay plots |
| `hi-simulate`    | PC      | fakes a full raw session (sky + noise + bandpass + RFI) so you can test everything before the telescope is done |
| `hi-stack`       | PC      | combines several processed sessions into one deeper RA map |
| `hi-check`       | both    | hardware bench test: ADC health, bias tee, bandpass, noise vs. radiometer limit |
| `hi-spectrum`    | PC      | analyze any single HI spectrum CSV (e.g. from the PICTOR web telescope): baseline, Gaussians, LSR, HI4PI overlay + Kelvin calibration |
| `hi-sunscan`     | PC      | solar pointing calibration: fits the Sun's drift through the beam from a daytime session |

## Quick start (no telescope needed yet)

```bash
uv sync
# simulate 24 h of drift scanning with a deliberately wrong pointing:
uv run hi-simulate --duration-h 24 --true-daz 5 --true-del -4
uv run hi-process data/raw/<session>
uv run hi-compare data/raw/<session> --strip hi4pi_cache/toy_strip.npz
```

`hi-compare` should report `daz ≈ +5, del ≈ -4` — it recovered the injected
pointing error from the data alone. The same procedure works on real data
against real HI4PI.

## Raspberry Pi setup

The **RTL-SDR Blog V4 requires the rtl-sdr-blog fork of the driver** — the
`rtl-sdr` package in Raspberry Pi OS is too old and misbehaves with the V4
(this, plus the Python 3.12 `distutils` removal, is why pyrtlsdr kept
failing; this pipeline uses neither — it talks to `rtl_tcp` over a socket
with pure Python).

```bash
sudo apt purge ^librtlsdr rtl-sdr
sudo apt install git cmake build-essential libusb-1.0-0-dev
git clone https://github.com/rtlsdrblog/rtl-sdr-blog
cd rtl-sdr-blog && mkdir build && cd build
cmake ../ -DINSTALL_UDEV_RULES=ON && make -j4
sudo make install && sudo cp ../rtl-sdr.rules /etc/udev/rules.d/ && sudo ldconfig
# verify: rtl_test  (should identify the V4 / R828D)
```

Then on the Pi:

```bash
git clone <this repo> && cd pipeline1
uv sync            # or: pip install -e .  (only numpy/pyyaml are needed to observe)
rtl_tcp -a 127.0.0.1 -p 1234 &
uv run hi-observe --config config.yaml           # Ctrl+C to stop
```

For unattended operation install the systemd units in `extras/`
(`systemctl enable --now rtl_tcp hi-observe`). The observer survives driver
hiccups (auto-reconnect) and writes self-contained session directories under
`paths.raw_dir` (set it to an SSD mount, e.g. `/mnt/ssd/hi-data`, in
`config.yaml`). Acquisition and processing are fully decoupled: sessions are
portable folders, so record on the Pi for days, then copy to any PC
(`rsync -a /mnt/ssd/hi-data/ laptop:.../data/raw/`) and run `hi-process`
there — as often as you like with different settings. What is stored per
15 s cycle is the averaged ON and OFF power spectrum at full 0.124 km/s
resolution (~32 KB): every downstream step (calibration, RFI, baselines,
velocity gridding) is recomputed from these on each `hi-process` run. Disk
budget: ~100 MB/day — a small SSD holds months.

If the LNA is powered externally keep `bias_tee: false` (the default here);
`hi-observe` and `hi-check` will then never put DC on the coax.

**Before first light, edit `config.yaml`:** site lat/lon, and your best guess
of the dish azimuth/elevation. A compass + inclinometer guess is fine —
`hi-compare` fixes it afterwards.

### Gain: the one knob that matters

The SAWbird already provides ~40 dB of gain. Start with `gain_db: 29.7`; run
`rtl_test` or watch `hi-observe`'s power printout — if spectra saturate
(power doesn't change when you disconnect the antenna, or ADC clips), lower
it. The HI line is detected in the *ratio* of ON/OFF spectra, so absolute
level only matters for staying in the ADC's linear range.

### Amplifier chain (SAWbird H1 vs. adding an LNA4ALL)

Friis: the first amplifier sets the system noise. The SAWbird H1
(NF ≈ 0.6 dB, ~40 dB gain) at the feed reduces everything behind it by a
factor 10⁴ — the RTL-SDR's own noise (~1000 K) contributes only ~0.1 K
through it. A second LNA therefore adds **no sensitivity**, and 60 dB of
total gain risks compressing the RTL-SDR front end, forcing its gain down.
Keep the chain minimal: `cantenna → SAWbird H1 → coax → RTL-SDR V4`.
Exceptions where the LNA4ALL earns its place:
* a very long/lossy coax run (≳ 15–20 dB loss) between SAWbird and SDR —
  insert the LNA4ALL *after* the long cable, right at the SDR;
* strong out-of-band RFI: insert the passive H-line filter (not the LNA)
  between SAWbird and SDR — its ~2 dB insertion loss is invisible after
  40 dB of gain, and `hi-check`'s quotient-noise figure tells you whether
  it actually helped.

## How the signal processing works

1. **Frequency switching** — every ~15 s the SDR hops between
   1420.1 MHz (ON: line in band) and 1423.2 MHz (OFF: line-free).
   `(ON − OFF) / smooth(OFF)` cancels the receiver bandpass and additive
   spurs (including the DC spike); multiplied by `T_sys` this gives antenna
   temperature. Without this step the HI line (~1% of system power) is
   invisible under the bandpass ripple.
2. **RFI flagging** — deviant cycles (broadband), persistent carriers
   (per-channel variance) and bursts (per-pixel outliers) are excised in
   *topocentric* channel space, where terrestrial RFI is stationary. The
   velocity window where galactic HI lives is protected from flagging.
3. **Doppler/LSR correction** — each cycle's frequency axis is shifted to the
   Local Standard of Rest (barycentric correction + 20 km/s solar apex
   motion, via astropy) and regridded to a common velocity grid; this also
   makes RFI smear out while the sky adds coherently.
4. **Baseline removal** — sigma-clipped polynomial per spectrum (asymmetric
   clipping so emission doesn't bias the fit).
5. **Products** — time×velocity waterfall, RA×velocity map, per-RA-bin:
   peak T_B, moments, **HI column density** (N_HI = 1.823·10¹⁸ ∫T_B dv),
   multi-Gaussian decomposition (spiral-arm components), and — where the
   scan crosses the galactic plane — **kinematic distances** (face-on
   Milky Way map) and a tangent-point **rotation curve** if the strip cuts
   the inner Galaxy.

### Sensitivity

Radiometer equation with T_sys ≈ 60 K, 1 km/s channels (4.7 kHz): a single
15 s cycle reaches σ ≈ 0.5 K — the galactic plane (30–100 K) is visible in
*one cycle*. A full night per 1° RA bin reaches ~0.05 K. Stack several days
of drifts at the same elevation to go deeper.

## HI4PI comparison & pointing calibration

HI4PI (Effelsberg + Parkes all-sky HI survey, 16′ resolution, in Kelvin) is
the perfect reference: your 80 cm dish sees the same sky through a ~17° beam.

```bash
uv run hi-fetch-hi4pi                  # downloads tiles for your dec strip (~5–9 GB once)
uv run hi-process data/raw/<session>
uv run hi-compare data/raw/<session>
```

`hi-fetch-hi4pi` grabs only the 20°×20° CDS tiles overlapping your strip and
reduces them to a ~50 MB cache (0.5°, 2.6 km/s — plenty under a 17° beam).
`hi-compare` then:

* smooths HI4PI to your beam, simulates what *your* telescope should have
  seen for a grid of azimuth/elevation offsets,
* cross-correlates each with your waterfall → the peak is your **true
  pointing** (a fixed az/el drift scan maps to a constant-declination track,
  so plane crossings + velocity structure pin it down to ≲1°, far below the
  beam width),
* the best-fit amplitude ratio gives your **effective T_sys**, i.e. an
  absolute Kelvin calibration for free,
* and writes side-by-side / overlay / spectra-comparison plots.

Update `pointing:` and `tsys_assumed_k:` in `config.yaml` with the fitted
values and re-run `hi-process` — now your maps are in true coordinates and
true Kelvin.

## Data formats

* raw session: `data/raw/<UTC>_<tag>/meta.json` + `chunk_*.npz`
  (`t_mid`, `on`, `off` = averaged |FFT|² spectra per cycle)
* products: `calibrated.npz` (t, v_LSR, T_B, ra, dec), `ra_map.npz`,
  `extraction.csv`, `pointing_fit.json`, PNG plots — all in
  `data/raw/<session>/products/`

## Practical tips

* Point the dish somewhere between south and the zenith-ish: from Lithuania
  the Milky Way (Cygnus/Cassiopeia region) drifts through elevations
  40–70° — rich HI and two galactic-plane crossings per day.
* Give it ≥ 24 h per session so `hi-compare` sees a full RA circle;
  the pointing fit gets much tighter with ≥ 1 plane crossing.
* Set `ppm` in `config.yaml` if your dongle has a known offset (V4s with
  TCXO are usually < 1 ppm — fine as is; 1 ppm = 1.4 kHz = 0.3 km/s).
* Check `rtl_biast`-style bias-tee behavior: `hi-observe` enables the
  bias tee via `rtl_tcp` on connect (`bias_tee: true`).
