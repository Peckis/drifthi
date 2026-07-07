# Roadmap — recommended additions, in rough priority order

What to build next, and why. Effort: S = an evening, M = a weekend, L = a project.

## 1. Nightly automation on the Pi (S)
The pipeline is manual-friendly but the telescope runs 24/7. Add:
* a systemd timer that restarts `hi-observe` at local noon daily → one tidy
  session per night, natural chunking for stacking;
* a cron job on the Pi: `hi-process --all && hi-stack --auto && hi-skymap`
  after each restart, then `rsync`/push products to the laptop (or commit
  small PNGs to a branch).
Result: wake up to processed maps without touching anything.

## 2. `hi-status` + failure notifications (S)
One command (and cron email/Telegram hook) reporting: is hi-observe alive,
minutes since last chunk, disk free, last night's cycle count, RFI %,
median P drift. The failure mode to defend against is silent: the service
dies at 2 am and you lose a week before noticing.

## 3. Continuum products (S–M, high science/effort ratio)
The broadband power per cycle (the `P` column) is already in the chunks —
currently only the sunscan uses it. Add a continuum branch to hi-process:
per-cycle total power → RA profile per session → stacked continuum drift
curves. Free science: Sun, the galactic-plane continuum ridge, and likely
Cas A (dec +58.8, ~1700 Jy — well placed for Lithuania) as a repeatable
calibrator transit. Also a great RFI/gain-stability diagnostic.

## 4. In-band frequency switching (M, +41% sensitivity)
Currently OFF (1423.2 MHz) contains no line — half the integration time
buys only bandpass calibration. Classic upgrade: shift by ~±0.7 MHz so the
line stays in band at BOTH tunings, then "fold" the two shifted copies.
Same bandpass cancellation, but all time is on-line: √2 better sensitivity.
Needs: config for the shift, a folding step in calibrate.py, care with the
DC-spike region overlapping the line at one tuning. Simulate first
(hi-simulate already models everything needed).

## 5. CI smoke test (S)
GitHub Actions: on push, run `hi-simulate --duration-h 0.5` →
`hi-process` → assert products exist and the recovered line amplitude is
sane. The end-to-end sim is already the de-facto test — make it automatic
so refactors can't silently break calibration.

## 6. Beam mapping from multiple sun transits (M)
One transit/day at a slightly different elevation offset each day for a
week → amplitudes + widths of each bump → a real 2-D beam profile (width,
ellipticity, pointing vs elevation). Improves hi-compare (true beam instead
of Gaussian) and quantifies the offset-feed geometry once and for all.

## 7. Observing campaigns (planning, not code)
* **Elevation ladder**: full sky ring every ~8° of elevation → your own
  all-sky HI map in ~2 months of unattended drift scanning (hi-skymap
  tracks progress).
* **Rotation curve**: low elevation (~20–30°) toward the south reaches
  dec +5..+20 → sightlines into the inner Galaxy (l ≈ 40–70°) →
  tangent-point rotation curve for R ≈ 5–8 kpc (extract.py already
  computes it when the data covers those longitudes).
* **HVC hunt**: deep stacks (many nights per ring) + extend the velocity
  grid (`v_min/v_max` to ±400) → high-velocity clouds at 1–3 K.

## 8. PICTOR email-verification loop (M)
The bot can't detect a swallowed slot. Poll the inbox via IMAP (Gmail app
password) for confirmation/CSV emails; resubmit missing slots at window
end; optionally auto-download the CSVs and run hi-spectrum on arrival.

## 9. Polyphase filterbank (L, polish)
Replace Hann+FFT with a 4-tap PFB in observe.py: much lower spectral
leakage → narrow RFI stays in one channel instead of splattering. Only
worth it if RFI at the site proves annoying; CPU cost ~2× (Pi 4 can take it).

## 10. Absolute calibration via Y-factor (M, hardware)
Point at ground/vegetation (~300 K) vs cold zenith sky (~10 K) at 1420 MHz;
the power ratio gives T_sys independently of HI4PI. Cross-check the
hi-compare gain fit; catches slow LNA degradation over seasons.

## Deliberately not recommended
* Bigger FFTs / finer channels — the 0.124 km/s raw channels already
  oversample any galactic HI line by ~50×.
* GPU processing — a full day reprocesses in under a minute on the laptop.
* pyrtlsdr — no.
