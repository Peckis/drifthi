# drifthi — the complete guide

How the pipeline works, module by module, and how to run it day to day.

---

## 0. The physics you're exploiting

**The line.** Neutral hydrogen emits at 1420.405751768 MHz (21 cm) from the
hyperfine spin-flip transition. The Milky Way is full of HI, so *every*
direction shows emission — the question is only its strength and shape. Gas
moving along the line of sight Doppler-shifts the line: 1 km/s ≈ 4.74 kHz.
Because the Galaxy rotates differentially, gas at different distances along
one sightline lands at different velocities — a single spectrum is a
one-dimensional slice through Galactic structure. Brightness: 1–100 K of
antenna temperature depending on direction; line widths tens of km/s.

**The drift scan.** Your dish is fixed at one azimuth/elevation. A direction
fixed to the Earth keeps a constant angle to the rotation axis, so the beam
traces a circle of **constant declination**, with right ascension advancing
at the sidereal rate. One sidereal day = one complete ring of sky. Change
the elevation between days to scan different declination rings. From
Lithuania (lat ≈ 55°), pointing due south at elevation *el* gives
dec ≈ *el* + 55° − 90°.

**The enemy: your own receiver.** The HI line is ~1% of total received
power, but the receiver bandpass (SAW filter + tuner + ADC response) ripples
by several percent across the band. In a raw spectrum the line is buried
under instrument shape. Everything in this pipeline exists to defeat that.

**The weapon: frequency switching.** The bandpass is fixed relative to the
tuner's local oscillator; the sky is fixed in absolute frequency. Retune the
LO and the bandpass moves with it while the sky doesn't. So we alternate:

* **ON**: centered 1420.1 MHz — line in band (at +0.31 MHz, away from the
  DC spike at band center)
* **OFF**: centered 1423.2 MHz — no line anywhere in band, still inside the
  SAWbird passband

and form the quotient (see §3). The result is the line spectrum divided by
system temperature, with the bandpass and additive spurs cancelled.

**Sensitivity (radiometer equation).** σ_T = T_sys / √(Δν · τ).
With T_sys ≈ 60 K on a 1 km/s channel (4.75 kHz): one 15 s cycle → σ ≈ 0.5 K
(galactic plane visible in a single cycle); one hour per RA bin → ≈ 0.06 K.
Noise integrates down as √time, forever, as long as the calibration is good —
which is what frequency switching guarantees.

---

## 1. The data flow

```
                     Raspberry Pi                                laptop
  ┌──────────┐   ┌───────────────┐   session folder    ┌──────────────────────┐
  │ rtl_tcp   │◄──┤ hi-observe    ├──► meta.json        │ hi-process           │
  │ (driver)  │   │ FFT+integrate │    chunk_00001.npz ─┼─► calibrated.npz     │
  └──────────┘   │ ON/OFF switch │    chunk_00002.npz  │   ra_map.npz         │
                  └───────────────┘    ...        rsync │   extraction.csv     │
                                                        │   *.png              │
  ┌────────────┐  one-time, ~5-9 GB                     └─────────┬────────────┘
  │ hi-fetch-  ├──► hi4pi_cache/strip_cache.npz (~50 MB)          │
  │ hi4pi      │                                        ┌─────────▼────────────┐
  └────────────┘                                        │ hi-compare           │
                                                        │ pointing + T_sys fit │
  hi-simulate ──► fake session folders (same format)    │ side-by-side plots   │
  hi-check    ──► hardware bench report                 └──────────────────────┘
```

Acquisition and processing are fully decoupled: a session folder is
self-describing (meta.json carries every parameter), so you can reprocess
years later with different settings.

---

## 2. Acquisition — `observe.py`, `sdr_tcp.py`, `session.py`

**`sdr_tcp.py`** is a ~80-line pure-Python client for the `rtl_tcp` protocol:
a TCP socket carrying 5-byte commands (set frequency, sample rate, gain,
bias tee...) one way and raw interleaved uint8 I/Q samples the other way.
This is why there's no pyrtlsdr and nothing to compile — any machine that
can run `rtl_tcp` (Pi, Windows, anything) can be observed from.

**`observe.py`** runs the eternal loop. Per cycle:

1. tune ON, discard `settle_s` of samples (PLL settling),
2. read exactly `t_on_s × sample_rate` samples; in blocks: convert uint8 →
   complex64, chop into `nfft`-sample pieces, apply a Hann window, FFT,
   accumulate |FFT|² — at the defaults that's ~8,790 spectra averaged into
   one 4096-channel ON spectrum,
3. same at the OFF tuning,
4. append (timestamp, ON, OFF, counts) to the current chunk; every
   `chunk_cycles` cycles, write a compressed `.npz` **atomically** (write
   temp file, rename) so a power cut never corrupts data.

Any socket/driver error triggers close-reconnect-retune after 5 s; Ctrl+C
(or systemd stop) flushes the partial chunk and finalizes `meta.json`. The
per-cycle log line `P_on/P_off` is a health signal: it should hover near
1.000; drifts mean gain instability, jumps mean RFI.

Why these defaults:
* `nfft: 4096` → 586 Hz = 0.124 km/s channels. HI features are ≥ few km/s
  wide, so this oversamples deliberately — you can always average down in
  processing, never up.
* `t_on_s: 15` — short enough that gain drift within a cycle is negligible,
  long enough that retune/settle overhead (~1 s) is small.
* `gain_db: 29.7` — with the SAWbird's ~40 dB in front, more SDR gain risks
  ADC clipping. Valid values are the discrete steps `rtl_test` prints.

**`session.py`** defines the on-disk format and `load_session()`, which
concatenates all chunks time-sorted. ~32 KB per cycle ≈ 100 MB/day.

---

## 3. Calibration — `calibrate.py`

Write the measured power in channel *ch* as

```
P_on(ch)  = G(ch)·[T_sys(ch) + T_line(ch)] + S(ch)
P_off(ch) = G(ch)· T_sys(ch)              + S(ch)
```

where G is the bandpass gain (identical at both tunings because it's
LO-relative) and S is additive junk (DC spike, ADC spurs — also identical).
Then:

```
P_on − P_off           G·T_line            T_line
────────────────  =  ────────────  ≈  ─────────────
 smooth(P_off)        ≈ G·T_sys           T_sys
```

* The **subtraction** kills additive terms *exactly* (raw OFF, not smoothed,
  is used in the numerator for this reason).
* The **division** kills the bandpass shape. The denominator is smoothed
  (median filter + Savitzky-Golay) so its noise isn't injected into every
  channel; it's floored at 0.1% of the spectrum median so the steep rolloff
  shoulders can't cause a divide-by-zero blowup.
* Multiplying by `tsys_assumed_k` converts to Kelvin. This is the one
  uncalibrated number in the system — and `hi-compare` measures it for you
  (§7), after which you update the config and reprocess.

`channel_mask()` discards the filter rolloff (`edge_frac`, 8% each side) and
the channels around the DC spike. `fit_baseline()` removes what calibration
couldn't: slow residual ripples (the bandpass isn't *perfectly* identical at
tunings 3.1 MHz apart). It fits a low-order polynomial with **asymmetric**
sigma-clipping — points far *above* the fit (emission!) are excluded
aggressively, points below only mildly — so the baseline rides the noise
floor without eating the line.

---

## 4. RFI excision — `rfi.py`

Terrestrial interference is narrowband and fixed in *frequency*; the sky
signal is fixed in *velocity*. Flagging therefore happens in channel space
**before** Doppler regridding, in three passes:

1. **Cycle flagging** — cycles whose broadband power (median over line-free
   channels) deviates > 5σ (MAD-based) are dropped whole: arcing thermostats,
   gain jumps.
2. **Channel flagging** — channels whose variance *over time* far exceeds
   the radiometer expectation are killed for the whole session: persistent
   carriers.
3. **Pixel flagging** — individual (cycle, channel) outliers vs. the
   per-channel median: bursty RFI.

The subtlety: the HI line itself is a legitimate ~100σ outlier. Channels in
the **protected window** (|v_topo| < 220 km/s, where galactic HI can occur)
are exempt from passes 2–3, so the pipeline can't flag its own signal. RFI
that lands inside the window survives flagging but gets diluted at the next
step: after LSR regridding it smears across velocity bins while the sky adds
coherently.

All thresholds are `processing.rfi_*` in the config. Flagged pixels become
NaN and are ignored (not interpolated over) by everything downstream.

---

## 5. Coordinates & velocities — `velocity.py`

* **Channel → velocity:** v = c(f₀ − f)/f₀ (radio convention).
* **Where you're pointing:** astropy AltAz→ICRS with your site and the
  configured az/el gives RA(t), Dec(t) for each cycle. (IERS auto-download
  is disabled so an offline Pi never stalls; the arcsecond-level error this
  causes is nothing under a 17° beam.)
* **Topocentric → LSR:** your measured velocities contain Earth's rotation
  (±0.3 km/s), Earth's orbit (±30 km/s, changes through the year!) and the
  Sun's drift relative to the local gas. Correction = astropy barycentric
  correction + 20 km/s projected toward the solar apex (18h03m +30°, the
  kinematic LSR definition). Every catalog (including HI4PI) uses v_LSR, so
  this is what makes your data comparable to anything.

`process.py` computes one scalar correction per cycle, shifts that cycle's
velocity axis, and linearly re-interpolates onto the common grid
(`v_min_kms … v_max_kms` step `dv_kms`). Then baselines (§3), then two data
products:

* **calibrated.npz** — the full time × velocity waterfall T(t, v) plus the
  RA/Dec track: the input for `hi-compare`.
* **ra_map.npz** — cycles averaged into `ra_bin_deg`-wide RA bins: the
  science map. Over multiple days at one elevation, bins just accumulate
  and the noise keeps dropping.

---

## 6. Science extraction — `extract.py`

Per RA bin (all written to `extraction.csv`):

* **noise σ** from emission-free velocities (|v| > 180 km/s), MAD-based;
* **emission mask**: channels > 3σ, slightly dilated;
* **moments**: mom0 = ∫T dv (K km/s), mom1 = intensity-weighted mean
  velocity, mom2 = velocity dispersion;
* **column density**: N_HI = 1.823×10¹⁸ · ∫T_B dv cm⁻² — exact if the gas
  is optically thin (good above |b| ≈ few degrees; an underestimate right
  in the plane);
* **Gaussian decomposition**: up to 4 components, found by iterative
  peak-subtract then refit jointly. Distinct components usually mean
  distinct structures along the sightline: local gas near 0 km/s, then
  spiral arms at increasingly negative velocities (in the outer Galaxy).

Globally, using galactic (l, b) computed for each bin:

* **Kinematic distances / face-on map.** Assume a flat rotation curve
  (V = 220 km/s, R₀ = 8.178 kpc). Differential rotation gives
  v_lsr = V₀ sin l (R₀/R − 1), invertible to the galactocentric radius R of
  any emission at velocity v. In the outer Galaxy (R > R₀) the line-of-sight
  distance is then unique: d = R₀cos l + √(R² − R₀²sin²l). Every bright
  (l, v) pixel becomes a point in the galactic plane → `faceon_map.png`, a
  top-down map of the Milky Way's gas *measured from your backyard*. (Your
  northern declination strips look mostly outward, which is exactly the
  regime without distance ambiguity.)
* **Rotation curve.** Where the strip cuts the inner Galaxy (0 < l < 90°),
  the highest-velocity gas on each sightline sits at the tangent point
  R = R₀ sin l, giving V(R) = v_terminal + V₀ sin l directly — the classic
  flat-rotation-curve / dark-matter measurement.

---

## 7. The HI4PI machinery — `hi4pi.py`, `compare.py`

**`hi-fetch-hi4pi`** downloads reference data. HI4PI (Effelsberg + Parkes,
2016) is the definitive all-sky HI survey: 16.2′ resolution, 1.29 km/s
channels, brightness temperature in Kelvin — i.e. exactly your observable,
just with a 4000× smaller beam. CDS distributes it as 20°×20° tiles
(`CAR_<row A–I><col 01–18>.fits`, ~252 MB each); the fetcher figures out
which rows your declination strip touches, downloads those tiles (once), and
bins them to a **strip cache**: 0.5° pixels, 2.58 km/s channels, ~50 MB.
Under a 17° beam nothing finer carries information.

**`hi-compare`** answers "what *should* my telescope have seen?" and fits
the difference:

1. Smooth the strip cache with a Gaussian of your `beam_fwhm_deg`
   (solid-angle weighted, RA-wrap safe) — HI4PI now looks the way your dish
   sees the sky.
2. For a trial pointing error (Δaz, Δel): recompute the drift track,
   interpolate model spectra along it → a model waterfall M(t, v) on the
   same grid as your observed O(t, v).
3. Score with Pearson correlation over all (t, v) pixels. Repeat over a
   coarse grid (±15° in 3° steps), then a fine grid (0.75° steps) around
   the peak.
4. At the best offset, the least-squares amplitude ratio
   gain = ΣOM/ΣM² calibrates your intensity scale:
   **T_sys,true = T_sys,assumed / gain** — an absolute calibration with no
   hot/cold load, courtesy of a 100 m telescope.

Why it works so well: a fixed-az/el scan has only two real unknowns (the
track's declination and its hour-angle offset — (Δaz, Δel) map onto these),
and galactic-plane crossings + velocity structure are extremely distinctive
along a 360° track. Validated in simulation: injected (+5.0°, −4.0°) error
recovered to (+5.25°, −3.75°), T_sys 60 K recovered as 61 K, with 8 h of
data. Outputs: `pointing_fit.json`, side-by-side waterfalls, contour
overlay, correlation maps, per-time spectra overlays.

**The feedback loop:** run `hi-compare` → put `az_fitted/el_fitted` into
`pointing:` and `tsys_implied_k` into `tsys_assumed_k` → rerun `hi-process`
→ maps are now in true coordinates and true Kelvin. The pointing offset is
a property of your mount setup, so it holds for later sessions at the same
dish setting.

---

## 8. Rehearsal & diagnostics — `simulate.py`, `check.py`

**`hi-simulate`** is a full forward model: takes a sky (HI4PI strip cache if
present, else a built-in toy Milky Way), pushes it through the same beam
smoothing and track logic, then applies everything real hardware does —
bandpass ripple + edge rolloff, slow gain drift, DC spike, radiometer noise
scaled to the actual integration time, a persistent carrier, random RFI
bursts — and writes a byte-identical session folder. Because it *shares the
sky code with the comparison* but the corruption code is independent of the
correction code, recovering an injected `--true-daz/--true-del` is a
genuine end-to-end test of the whole chain.

**`hi-check`** is the hardware bench test (no dish needed): ADC health
(mean/std/clipping → is the gain right?), optional bias-tee power test
(disabled while your LNA is externally powered), ON/OFF bandpass plots, and
the killer metric: **quotient noise vs. the radiometer prediction**
√(2/(B_eff·τ)). At the limit → the whole chain works; far above → RFI or
instability. Use it to A/B test hardware changes (extra filter, different
gain, different location) with numbers instead of vibes.

---

## 9. Cookbook

**Bench (now, no dish):**
```
tools\rtl-sdr\x64\rtl_tcp.exe -a 127.0.0.1 -p 1234     # terminal 1
uv run hi-check                                         # terminal 2
```
Healthy = no clipping, quotient near radiometer limit, ADC std jumps up
when the SAWbird powers on.

**Rehearse the full loop (no hardware at all):**
```
uv run hi-simulate --duration-h 24 --true-daz 5 --true-del -4
uv run hi-process  data/raw/<session>
uv run hi-compare  data/raw/<session> --strip hi4pi_cache/toy_strip.npz
```

**Electronics first light (cantenna out a window, still no dish):** run
`hi-observe` overnight; the galactic plane should appear in the waterfall.

**Real operations:**
```
# once, on the laptop: reference data for your dec strip
uv run hi-fetch-hi4pi

# on the Pi (config: paths.raw_dir -> SSD, bias_tee: false):
systemctl enable --now rtl_tcp hi-observe        # units in extras/

# whenever you feel like science, on the laptop:
rsync -a pi:/mnt/ssd/hi-data/ data/raw/
uv run hi-process data/raw/<session>
uv run hi-compare data/raw/<session>
# first time: copy fitted az/el + T_sys into config.yaml, rerun hi-process
```

**Survey mode:** leave the dish for ≥ 1 full day per elevation, then change
elevation by ~half a beam (8–10°) and repeat — each setting adds a ring to
an all-sky HI map.

### First-light checklist

1. `config.yaml`: set `site` lat/lon and your best-guess `pointing` az/el.
2. Bench check with everything connected as it will observe:
   `rtl_tcp -a 127.0.0.1 -p 1234 -b 4` + `hi-check` (LNA powered
   externally, `bias_tee: false`). Want: no clipping, quotient noise near
   the radiometer limit, ADC std well above the bare-dongle value.
3. **Daytime, dish in its final position:** run
   `hi-observe --tag sunscan` from ~3 h before to ~3 h after solar noon,
   then `hi-sunscan data/raw/<session>_sunscan`. It fits the Sun's drift
   through the beam and reports the transit time, your **measured beam
   FWHM**, and the along-path pointing correction (a single transit cannot
   constrain the cross-path component -- that only lowers the bump height;
   `hi-compare` fixes it later, or tilt the dish ~5 deg and repeat).
   Sessions tagged `sunscan` are automatically excluded from `--all`/
   `--auto` science processing. Note: this needs the Sun's declination
   (+23 in July) within ~a beam of the dish declination.
4. Update `pointing` (and `beam_fwhm_deg`) from the sunscan, then observe
   all night, and every morning run the two-liner in the next section.
5. After the first full night: `hi-fetch-hi4pi` + `hi-compare` to nail the
   remaining pointing component and T_sys; put both into the config and
   reprocess. Done -- from then on it's just data accumulation.

### The daily two-liner (plus one for fun)

```bash
uv run hi-process --all     # processes only sessions with new data
uv run hi-stack --auto      # groups sessions by declination, restacks each
uv run hi-skymap            # all-sky Aitoff map: every track's beam smear,
                            # the Milky Way band, and % of sky / MW covered
```

`--all` skips up-to-date sessions (`--force` overrides); `--auto` groups
sessions whose declinations agree to `--dec-tol` (default 2 deg) and writes
one deepening stack per dish setting under `data/stacks/dec<+XX>/`. Every
session and stack gets a `coverage.png` -- the beam swept along the drift
track (an elongated band, not a circle: the sky moves during observing),
drawn over the HI4PI column-density map when a strip cache exists.

### Running as a service on the Pi

```bash
sudo cp extras/rtl_tcp.service extras/hi-observe.service /etc/systemd/system/
# edit hi-observe.service: User=, WorkingDirectory=, ExecStart= paths
sudo systemctl daemon-reload
sudo systemctl enable --now rtl_tcp hi-observe   # start now + at every boot

systemctl status hi-observe          # is it running? last log lines
journalctl -u hi-observe -f          # live log (cycle counter, P_on/P_off)
sudo systemctl stop hi-observe       # clean stop
sudo systemctl restart hi-observe    # e.g. after editing config.yaml
```

`stop` sends SIGTERM; the observer treats it exactly like Ctrl+C — it
flushes the buffered cycles into a final chunk and finalizes `meta.json`.
Chunks are also written atomically every ~10 min, so even a power cut only
costs the unwritten buffer, never corrupts a session. Each (re)start opens
a *new* session folder; that's normal — `hi-stack` merges them later.

### Getting data to the laptop

`rsync` copies only what's new: it compares files on both ends and skips
chunks it already transferred, so syncing a week of data after a day of new
observing transfers one day. Safe to run *while observing* (chunks are
atomic; the live session simply gains files on the next sync). From Linux/
WSL/Git-Bash with rsync installed:

```bash
rsync -av pi@raspberrypi.local:/mnt/ssd/hi-data/ data/raw/
```

(`-a` = preserve everything, recurse; `-v` = show files; trailing slashes =
"contents into contents".) On plain Windows, OpenSSH's `scp -r` works but
recopies everything; WinSCP's "synchronize" gives rsync-like behavior with
a GUI.

### Partial processing and stacking

```bash
uv run hi-process data/raw/<session> --last-h 2      # quick look: newest 2 h
uv run hi-process data/raw/<session> --start-h 3 --stop-h 9
uv run hi-stack  data/raw/night1 data/raw/night2 --out data/stacks/week1
```

`hi-process` is stateless: it always recomputes everything from the raw
chunks and overwrites `products/`. Re-running after more data arrived just
means running it again — no incremental bookkeeping to get wrong. To
combine *multiple sessions* (several nights, or restarts) use `hi-stack`:
it averages their RA maps weighted by cycle counts (noise drops ~sqrt(N)
where they overlap, gaps fill where they don't) and re-runs the science
extraction on the deep map. Only stack sessions with the same dish
pointing.

### Analyzing spectra from other telescopes (PICTOR etc.)

`hi-spectrum` runs the single-spectrum half of the pipeline on any CSV of
frequency vs. power — e.g. the file PICTOR emails you:

```
uv run hi-spectrum observation.csv --pictor --l 120 --b 5 \
    --time 2026-07-10T21:30:00 --fetch-hi4pi
```

It auto-detects the CSV layout and frequency units, applies the LSR
correction for the site/time/pointing, fits and removes a baseline, finds
Gaussian components, and — if HI4PI data covers the pointing
(`--fetch-hi4pi` downloads the 1–4 needed tiles) — overlays the
beam-matched model and fits the amplitude scale, converting the arbitrary
power units to Kelvin so N_HI comes out in cm⁻². What it can *not* do for
foreign data is bandpass calibration or RFI excision — those need the
ON/OFF cycles and time series that only exist in our own raw sessions.

## 10. Troubleshooting

| symptom | likely cause | fix |
|---|---|---|
| `hi-check`: ADC rail % > 0 | too much gain | lower `sdr.gain_db` one step |
| ADC std < 3 | LNA unpowered / cable | check SAWbird LED, connections |
| quotient noise ≫ radiometer limit | RFI, USB dropouts, gain instability | try H-line filter, shorter USB, other location |
| waterfall shows vertical stripes | persistent RFI carrier | raise `rfi_chan_sigma`? usually flagged automatically; check bandpass plot |
| waterfall horizontal stripes | bad cycles slipping through | lower `rfi_cycle_sigma` |
| broad curved residuals in spectra | baseline order too low | `baseline_order: 5` |
| line clipped at map edge | velocity grid too narrow | widen `v_min/v_max_kms` |
| hi-compare r < 0.3 | wrong site coords, dec strip not in cache, or too-short session | check lat/lon sign, refetch strip, observe ≥ 12 h |
| hi-compare peak at search edge | true error > search range | raise `compare.search_deg` |
