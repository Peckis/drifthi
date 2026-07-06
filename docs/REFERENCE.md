# drifthi reference manual

Exhaustive documentation of every command, option, config key, and output
file. For concepts and physics see [GUIDE.md](GUIDE.md); this file is the
boring, complete "what does this knob do" list.

All commands are installed by `uv sync` and run as `uv run <command>` (or
directly from `.venv/bin/` / `.venv\Scripts\`). Every command accepts
`--config PATH` (default `config.yaml` in the current directory) and
`-h/--help`. All timestamps everywhere are **UTC**. Ctrl+C is always safe:
`hi-observe` flushes its buffer and finalizes the session before exiting;
all other tools are pure analysis and can be re-run at will.

---

## 1. config.yaml — every key

### `paths`
| key | default | meaning |
|---|---|---|
| `raw_dir` | `data/raw` | where `hi-observe`/`hi-simulate` create session directories, and where `--all`/`--auto` modes look. Point at an SSD mount on the Pi. |

### `site`
| key | meaning |
|---|---|
| `name` | free text, recorded in session metadata |
| `lat_deg` | geographic latitude, degrees north |
| `lon_deg` | geographic longitude, degrees **east positive** |
| `height_m` | elevation above sea level (low impact; rough is fine) |

### `pointing`
| key | meaning |
|---|---|
| `az_deg` | dish azimuth, degrees clockwise from north (180 = due south) |
| `el_deg` | dish elevation above horizon |
| `beam_fwhm_deg` | beam full width at half maximum. Initial estimate ~1.22·λ/D; **hi-sunscan measures the real value** — update it here. Used by beam smoothing, coverage plots, sunscan geometry. |

The resulting drift declination is `dec ≈ el + lat − 90` when az=180;
`hi-observe` computes it exactly per cycle. `hi-compare` fits the true
az/el — feed its result back here.

### `sdr`
| key | default | meaning |
|---|---|---|
| `host`, `port` | 127.0.0.1, 1234 | where rtl_tcp listens. If nothing is listening on localhost, `hi-observe`/`hi-check` start `rtl_tcp -b 4` themselves (and stop it on exit). |
| `sample_rate_hz` | 2.4e6 | SDR sample rate = usable bandwidth. 2.4 MS/s is the V4's comfortable maximum. |
| `freq_on_hz` | 1420.1e6 | ON tuning: HI line lands +0.31 MHz from center, away from the DC spike. |
| `freq_off_hz` | 1423.2e6 | OFF tuning: line-free, still inside the SAWbird passband. Must differ from ON by ≲ a few MHz (bandpass similarity) and ≳ 1.5 MHz (no line overlap). |
| `gain_db` | 29.7 | tuner gain. Only the discrete steps printed by `rtl_test` exist; invalid values are rounded by the driver. Judge with `hi-check`: raise while ADC std < ~3, lower on any clipping. |
| `ppm` | 0 | frequency correction. V4's TCXO is <1 ppm; 1 ppm = 0.3 km/s. |
| `bias_tee` | true | power the LNA through the coax. **false** if powered externally — then no DC ever goes up the cable, and `hi-check` skips its bias test. |
| `nfft` | 4096 | FFT length → channel width fs/nfft (586 Hz = 0.124 km/s). Larger = finer channels, more CPU. |
| `t_on_s`, `t_off_s` | 15, 15 | integration per tuning per cycle. Shorter = more retune overhead; longer = more gain drift within a cycle. |
| `settle_s` | 0.5 | data discarded after each retune. Must exceed PLL settling *and* the rtl_tcp ring buffer (hence `-b 4`). Don't reduce. |
| `chunk_cycles` | 20 | cycles per output file (~10 min). Only affects file granularity / max data at risk on power loss. |

### `processing`
| key | default | meaning |
|---|---|---|
| `tsys_assumed_k` | 60 | system temperature used to scale the quotient into Kelvin. A placeholder until `hi-compare` measures the real value — then update and reprocess. |
| `edge_frac` | 0.08 | fraction of channels dropped at each band edge (filter rolloff) |
| `dc_halfwidth_bins` | 3 | channels blanked around the center (DC spike residue) |
| `v_min_kms`, `v_max_kms`, `dv_kms` | −250, +250, 1.0 | the common LSR velocity grid all spectra are interpolated onto |
| `baseline_order` | 3 | polynomial order of the per-spectrum baseline. Raise to 4–5 if broad curved residuals remain. |
| `rfi_cycle_sigma` | 5 | MAD threshold for dropping whole cycles (broadband events) |
| `rfi_chan_sigma` | 6 | threshold for killing persistently noisy channels |
| `rfi_pixel_sigma` | 6 | threshold for individual (cycle, channel) outliers |
| `ra_bin_deg` | 1.0 | RA bin width of the science map |
| `time_bin_s` | 240 | time binning used by `hi-compare` |

### `hi4pi`
| key | default | meaning |
|---|---|---|
| `cache_dir` | `hi4pi_cache` | where survey tiles and strip caches live |
| `base_url` | CDS EQ2000/CAR | tile download location |
| `v_max_kms` | 350 | velocity crop of the strip cache |
| `sky_bin_deg` | 0.5 | strip cache spatial resolution |
| `v_bin_ch` | 2 | survey channels averaged per cache channel (2 → 2.58 km/s) |

### `compare`
| key | default | meaning |
|---|---|---|
| `search_deg` | 15 | ± range of the pointing-offset search |
| `coarse_step_deg` | 3 | first-pass grid step |
| `fine_step_deg` | 0.75 | refinement grid step |

---

## 2. Commands

### hi-observe — record a drift-scan session (runs on the Pi)

```
hi-observe [--config C] [--out DIR] [--tag TAG] [--duration-h H] [--no-bias-tee]
```
| option | default | meaning |
|---|---|---|
| `--out DIR` | `paths.raw_dir` | root under which the session directory is created |
| `--tag TAG` | `obs` | session name suffix → `<UTC>_<tag>`. Use `--tag sunscan` for solar sessions: those are excluded from `--all`/`--auto` science processing. |
| `--duration-h H` | 0 = forever | stop automatically after H hours |
| `--no-bias-tee` | | force the tee off for this run regardless of config |

Behavior: creates the session dir immediately; writes a chunk every
`chunk_cycles` cycles (atomic: temp file + rename); auto-reconnects on any
driver/socket error; starts rtl_tcp itself if none is running on localhost
(and stops it on exit). Stop with Ctrl+C or `systemctl stop hi-observe` —
both flush the buffer and finalize `meta.json`.

Log line anatomy:
`cycle 42 12:34:56Z P_on/P_off=1.0002 P=+0.15 dB (31.5s)`
— `P_on/P_off` should hover near 1 (jumps = RFI, drift = gain instability);
`P` is broadband power relative to the session start: the Sun shows up here
live as a smooth rise to ~+3 dB and back.

### hi-check — hardware bench test

```
hi-check [--config C] [--host H] [--port P] [--out DIR] [--bias-test] [--skip-bias-test]
```
| option | meaning |
|---|---|
| `--host/--port` | override `sdr.host`/`sdr.port` |
| `--out DIR` | report root (default `checks/`, one timestamped subdir per run) |
| `--bias-test` | force the bias off/on test even if `sdr.bias_tee` is false |
| `--skip-bias-test` | never toggle the tee |

Sequence: connect (auto-starting rtl_tcp if needed) → enable tee per
config → bias off/on power test (expect several dB if the LNA is
coax-powered) → ADC statistics in the operating state (mean ~127.5; std ≳3;
rails ~0%) → 4 s ON and OFF spectra → quotient noise vs the radiometer
prediction (~1× = perfect). Outputs `bench_spectra.png` + `report.json`.

### hi-sunscan — solar pointing calibration

```
hi-sunscan SESSION [--config C]          # analyze a daytime session
hi-sunscan --predict                     # when does the Sun cross the beam?
hi-sunscan --sun-at 2026-07-06T14:00:00  # where is the Sun at this UTC time?
```
`--predict` scans the next 48 h for the closest approach to the configured
pointing and rates the expected bump ("STRONG/detectable/too far"). Record
from ~2.5 h before to ~2.5 h after that time (minimum ~4 h total so the fit
has baseline on both sides; there is no maximum). `--sun-at` supports the
"calibrate at a convenient hour" trick: point the dish where the Sun will
be at that time (align by the feed shadow) and keep that pointing for the
night. The session analysis finds the bump automatically (no time input),
and reports: fitted transit time, **measured beam FWHM**, the along-path
pointing correction, and corrected az/el. One transit cannot constrain the
cross-path component (it only lowers the bump amplitude) — `hi-compare`
fixes that, or repeat after tilting the dish ~5° in elevation. Outputs
`products/sunscan.png` + `sunscan.json`.

### hi-process — raw session → calibrated products

```
hi-process SESSION [SESSION ...] [--all] [--force]
           [--start-h H] [--stop-h H] [--last-h H] [--config C]
```
| option | meaning |
|---|---|
| `--all` | process every session under `paths.raw_dir` that has chunks newer than its products (sunscan-tagged sessions skipped) |
| `--force` | with `--all`: reprocess even up-to-date sessions (use after changing config) |
| `--start-h/--stop-h H` | restrict to a time window, hours from session start |
| `--last-h H` | only the newest H hours (quick look at a running session — safe: chunks are atomic) |

Stateless: always recomputes everything from raw chunks and overwrites
`products/`. Pipeline: quotient calibration → RFI flagging → LSR
regridding → baselines → time & RA maps → extraction → plots.

### hi-stack — combine sessions into deeper maps

```
hi-stack SESSION [SESSION ...] [--out DIR]
hi-stack --auto [--dec-tol D]
```
| option | default | meaning |
|---|---|---|
| `--auto` | | discover all processed sessions, group by declination, write one stack per dish setting to `data/stacks/dec<+XX>/` |
| `--dec-tol D` | 2.0 | max dec difference for two sessions to share a group |
| `--out DIR` | timestamped | explicit output dir (non-auto mode) |

Cycle-weighted average per (RA bin, channel); overlaps get deeper (√N),
gaps fill in. Warns if declinations within one stack span > 3°. Requires
all inputs processed with the same velocity/RA grids (reprocess with
`--force` after grid changes).

### hi-skymap — all-sky coverage figure

```
hi-skymap [SESSION ...] [--out FILE] [--band-b B]
```
| option | default | meaning |
|---|---|---|
| positional | all processed sessions | restrict to specific sessions |
| `--out` | `data/skymap.png` | output image |
| `--band-b B` | 10 | \|b\| defining the "Milky Way band" for the statistics |

Aitoff projection, exact spherical beam smears, Milky Way band shaded,
galactic plane dashed, Galactic center starred; solid-angle percentages of
sky and band in the title.

### hi-fetch-hi4pi — download the reference survey

```
hi-fetch-hi4pi [--dec-min D] [--dec-max D] [--margin-deg M] [--delete-tiles]
```
| option | default | meaning |
|---|---|---|
| `--dec-min/--dec-max` | auto | explicit strip; default = pointing dec ± margin |
| `--margin-deg M` | beam + search range | half-width of the auto strip |
| `--delete-tiles` | keep | delete each 252 MB tile after binning (saves disk, costs re-download later) |

Downloads only the CDS tiles overlapping the strip (resumable; tiles are
cached, so re-runs skip finished files) and reduces them to
`hi4pi_cache/strip_cache.npz`. One-time per declination range.

### hi-compare — fit pointing + T_sys against HI4PI

```
hi-compare SESSION [--strip FILE] [--config C]
```
`--strip` selects a cache file (default `<cache_dir>/strip_cache.npz`).
Needs ≥ several hours of processed data; a full 24 h with a galactic-plane
crossing gives the tightest fit. Outputs the fitted az/el (put into
`pointing:`), implied T_sys (put into `tsys_assumed_k`, then reprocess),
`pointing_fit.json`, and the side-by-side / overlay / correlation plots.

### hi-simulate — synthesize a raw session

```
hi-simulate [--duration-h H] [--tag T] [--start-unix U] [--true-daz D]
            [--true-del D] [--strip FILE] [--rfi-prob P] [--seed N] [--out DIR]
```
`--true-daz/--true-del` inject a pointing error the pipeline should recover
(end-to-end test). Sky: `--strip` file if given, else the HI4PI cache if
present, else a built-in toy Milky Way. `--rfi-prob` is the per-cycle
probability of an RFI burst.

### hi-spectrum — analyze any single HI spectrum CSV

```
hi-spectrum FILE [--pictor] [--ra R --dec D | --l L --b B] [--time UTC]
            [--site-lat X --site-lon Y] [--beam-fwhm F] [--strip FILE]
            [--fetch-hi4pi] [--baseline-order N] [--edge-frac F]
            [--freq-col N] [--power-col N] [--center-mhz M --bandwidth-mhz B]
            [--out DIR]
```
| option | meaning |
|---|---|
| `--pictor` | use the PICTOR site (Athens) and 10° beam |
| `--ra/--dec` or `--l/--b` | pointing (equatorial or galactic) |
| `--time UTC` | observation time — needed for the LSR correction |
| `--fetch-hi4pi` | download the 1–4 tiles around the pointing for the model overlay + Kelvin calibration |
| `--edge-frac` | crop before baseline fitting (default 0.05; raise to 0.10–0.15 for steep-edged spectra) |
| `--freq-col/--power-col` | force CSV columns (default: first numeric = frequency, last = power) |
| `--center-mhz/--bandwidth-mhz` | build the frequency axis for single-column files |

### hi-pictor-bot — auto-observe with the PICTOR web telescope

```
hi-pictor-bot --email ADDR [--target cygnus|anticenter|both] [--window-h W]
              [--gap-s G] [--nights N] [--name PREFIX] [--dry-run]
```
| option | default | meaning |
|---|---|---|
| `--email` | required | where PICTOR sends CSVs (`raw_data=1` is always set) |
| `--target` | cygnus | which galactic-plane crossing of PICTOR's zenith to observe |
| `--window-h` | 1.0 | observing window centered on the crossing |
| `--gap-s` | 240 | pause between the 600 s requests (shorter risks the site silently dropping one) |
| `--nights` | 1 | repeat for N nights (the process sleeps in between — run under nohup/tmux) |
| `--dry-run` | | print schedule and would-be POSTs, submit nothing |

---

## 3. Data formats

**Raw session** `data/raw/<UTC>_<tag>/`:
* `meta.json` — every acquisition parameter + site + pointing snapshot
* `chunk_NNNNNN.npz` — `t_mid` (N,) unix; `on`, `off` (N, nfft) mean |FFT|²
  spectra (fftshifted, absolute power, arbitrary units); `n_on`, `n_off`
  (N,) FFT counts. ~32 KB/cycle, ~100 MB/day.

**Products** `<session>/products/` (all regenerated by every `hi-process` run):
* `calibrated.npz` — `t`, `v` (LSR grid), `T` (t×v Kelvin waterfall), `ra`,
  `dec`, `vcorr`, `meta` (JSON string)
* `ra_map.npz` — `ra`, `dec`, `l`, `b`, `v`, `W` (RA×v map), `counts`
* `extraction.csv` — per RA bin: coordinates, noise, peak, moments, N_HI,
  velocity extent, up to 4 Gaussian components (amp/v0/sigma each)
* `waterfall_time.png`, `waterfall_ra.png`, `avg_spectrum.png`,
  `column_density.png`, `coverage.png`; `faceon_map.png` and
  `rotation_curve.{png,csv}` when sky coverage allows
* from hi-compare: `pointing_fit.json`, `hi4pi_side_by_side.png`,
  `hi4pi_overlay.png`, `hi4pi_pointing_fit.png`, `hi4pi_spectra_overlay.png`
* from hi-sunscan: `sunscan.json`, `sunscan.png`

**Stacks** `data/stacks/<name>/`: same `ra_map.npz` + extraction + plots as
a session, with `counts` summed over all inputs.

**HI4PI strip cache** `hi4pi_cache/*.npz`: `ra`, `dec`, `v` (LSR km/s),
`cube` (ra×dec×v, Kelvin, 0.5° / 2.58 km/s).
