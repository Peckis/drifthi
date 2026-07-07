# The excruciatingly in-depth guide

A lecture in four parts: the physics from zero, what the software actually
does when you type a command, how to inspect and improve the physical
build, and the three planned upgrades explained properly. No prior radio
astronomy assumed. Read with coffee.

---

# Part I — Radio astronomy from zero

## 1.1 What the telescope actually measures: noise

Forget images. A radio telescope of this kind measures exactly one thing:
**how much random electrical noise comes out of the antenna, as a function
of frequency and time.** That's it. Every product this pipeline makes is a
rearrangement of that.

The universe at 1420 MHz doesn't send signals; it glows. Hydrogen gas,
the ground, trees, your receiver's own transistors — everything warm emits
random electromagnetic noise. Radio astronomers measure noise power in
temperature units, Kelvin, because of a beautiful fact: a resistor at
temperature T produces noise power P = k·T·B (Boltzmann constant times
bandwidth). So "this channel shows 5 K" means "the noise power here is
what a 5 K resistor would make." It's a currency, and everything trades
in it:

* **T_sys (system temperature)**: all the noise you *don't* want, added
  up — the LNA's transistors (~35 K for a SAWbird), losses in cables and
  connectors before the LNA (big! see Part III), ground radiation leaking
  into the feed (~tens of K), the atmosphere (~few K), the cosmic
  microwave background (2.7 K). For your build, expect T_sys ≈ 60–150 K.
* **T_A (antenna temperature)**: the signal — what the sky adds. Galactic
  hydrogen adds between ~1 K (looking out of the plane) and ~100 K
  (looking along the plane).

So the game is: detect a 1–100 K signal sitting on a 100 K floor of your
own noise. Sounds hopeless — it isn't, because of averaging (§1.4).

## 1.2 The 21 cm line, and why the spectrum has a shape

A neutral hydrogen atom is a proton and an electron; each has spin. When
the spins flip from parallel to antiparallel, the atom emits a photon at
exactly **1420.405751768 MHz**. A given atom does this about once per 10
million years — but a telescope beam through the Milky Way contains ~10^60
atoms, so the glow is steady and bright.

If all the gas were stationary you'd see one narrow spike at 1420.406.
Instead you see bumps and shoulders spread over ±1 MHz, because the gas
**moves**, and motion Doppler-shifts the line: 1 km/s of line-of-sight
velocity = 4.74 kHz of frequency shift. The Galaxy rotates differentially
(inner parts orbit faster), so gas at different distances along your line
of sight arrives at different velocities → different frequencies. **The
shape of the spectrum is a map of gas along the sightline.** That double
peak you saw in the PICTOR Cygnus data — nearby Local Arm gas at ~0 km/s
and Outer Arm gas at −75 km/s, 38,000 light-years away — is exactly this.

One habit to build: the pipeline always converts frequency to velocity
(v_LSR, "Local Standard of Rest" — corrected for Earth's spin, Earth's
orbit, and the Sun's drift, so measurements from any day and any telescope
line up). Negative velocity = gas approaching. Higher frequency = negative
velocity. You'll stop thinking in MHz within a week.

## 1.3 The antenna: what your dish + can actually do

An antenna doesn't amplify; it *selects direction*. Your 80 cm dish
focuses radiation from a patch of sky ~17° across (the **beam**) onto the
can (the **feed**), whose probe converts the electromagnetic field into a
voltage on the coax. Two facts about this that matter for your anxiety:

**Fact 1 — For extended emission, dish size does not set signal strength.**
If the emission fills the beam (galactic HI always does), the antenna
temperature equals the sky's brightness temperature, whatever the dish
diameter. A 100 m dish looking at the galactic plane and your 80 cm dish
looking at the galactic plane both register ~the same Kelvin. The big dish
sees a *smaller patch* (resolution), not a *stronger signal*. You are not
competing on sensitivity with Effelsberg; you're only competing on beam
size, and you already know your beam is ~17°.

**Fact 2 — Build imperfections cost efficiency, not detection.** A
misfocused feed, an imperfect can, a dented dish — these scatter some of
the power that should reach the probe, and let some ground noise sneak
in. The combined effect is a lower "beam efficiency" η: your measured line
is η × the true brightness temperature, and T_sys is somewhat higher.
Suppose your build is bad enough to cost a full factor of 2 (that's a
*very* pessimistic 3 dB): the galactic plane at 80 K then reads 40 K on
a 150 K system. Your noise per 15 s cycle is 150/√(4750×15) ≈ 0.6 K.
That's still a 70σ detection *per cycle*. Detection is not in question;
only the map depth per night is, and even that just means "stack more
nights." The calibration against HI4PI (hi-compare's gain fit) absorbs η
automatically — your Kelvin scale ends up right *even though* the build
is imperfect. This is why the pipeline was designed around HI4PI.

## 1.4 The radiometer equation: why averaging works

Noise is random, so it averages down. Measure a bandwidth Δν for a time τ
and you've collected N = Δν·τ independent samples of the noise; the
average of N random things fluctuates √N less than one of them:

    σ_T = T_sys / √(Δν · τ)

Your numbers: Δν = 4750 Hz (a 1 km/s channel), T_sys = say 100 K.
* one 15 s cycle: σ ≈ 0.37 K → the plane (30–100 K) is obvious instantly
* one hour: σ ≈ 0.048 K
* a whole night per 1° RA bin (~2 h effective): σ ≈ 0.035 K → high-latitude
  wisps at 1–3 K are easy

The catch — and this is the single most important idea in the whole
project — is that √N averaging only removes *random* error. Any
*systematic* error (a ripple in your electronics' frequency response)
does not average away, ever. It just gets measured more and more
precisely. Which brings us to:

## 1.5 Frequency switching, explained three ways

**The problem.** What the SDR records in channel f is
P(f) = G(f)·[T_sys + T_line(f)]: the sky *times the gain curve* G(f) of
the whole chain (filter ripple, tuner response, reflections in the coax).
G(f) wiggles by several percent across the band. The line is also a few
percent. Looking at raw P(f), you cannot tell wiggle from hydrogen — and
no amount of integration helps, because the wiggle is systematic (§1.4).
You measured this yourself with hi-check: 1.1–1.3 dB of ripple, i.e. the
instrument's shape is ~30×, the strongest possible line.

**Way 1 — the kitchen scale.** To weigh flour you put the empty bowl on
the scale and press "tare," then measure. The OFF spectrum is the tare:
the same instrument, measured with nothing (no line) in it. Division by
the OFF removes the instrument's shape. And because your scale *drifts*
(gain changes with temperature — you watched PICTOR's drift 9% in ten
minutes), you must re-tare constantly: hence switching every 15 s rather
than once per night.

**Way 2 — the algebra.** The trick works because of *where* the
bandpass lives. G(f) is created inside the receiver, *after* the mixer
that shifts the sky down by the local oscillator frequency. So G is glued
to the LO: retune the LO by 3 MHz and G stays put in channel space while
the sky slides by 3 MHz. Measure with the line in band (ON) and with the
LO shifted so no line is in band (OFF):

    P_on(ch)  = G(ch)·[T_sys + T_line(ch)] + spurs(ch)
    P_off(ch) = G(ch)· T_sys              + spurs(ch)

    (P_on − P_off) / smooth(P_off)  =  T_line / T_sys

Subtraction kills additive junk (DC spike, spurs — identical in both).
Division kills G exactly. What remains is the line in units of T_sys —
multiply by T_sys (calibrated via HI4PI) and you have Kelvin.

**Way 3 — what you'd see.** Plot a raw ON spectrum: a big lumpy hump,
line invisible. Plot ON and OFF on top of each other: two big lumpy humps
that look identical. Plot their difference: a flat line at zero with a
little bump sticking up at 1420.4. That bump is the Galaxy. (The
hi-check bench plot's two panels are literally this demonstration.)

**The cost.** Half your time goes to OFF, and the subtraction adds the
OFF's noise: net factor ~2 in observing time versus a hypothetical
perfect receiver. Perfect receivers don't exist; every professional
single-dish HI survey — including HI4PI itself — pays this same tax in
one form or another (they call it frequency switching, position
switching, or Dicke switching). Upgrade #2 in Part IV reduces the tax.

---

# Part II — What actually happens when you type a command

## 2.1 `hi-observe` — the recorder

What runs, step by step:

1. **Connects to rtl_tcp** (starting it if needed). rtl_tcp is a tiny
   server from the driver package that owns the USB dongle and streams
   raw samples over a local socket. We use it (rather than USB directly)
   because its protocol is 5-byte commands + a byte stream — nothing that
   can break.
2. **Configures the dongle**: sample rate 2.4 MHz (this *is* your
   bandwidth — the SDR delivers a 2.4 MHz-wide slice of spectrum centered
   wherever it's tuned), gain, bias tee.
3. **The cycle loop**, forever:
   a. tune to 1420.1 MHz (ON); throw away 0.5 s of samples (the tuner's
      synthesizer needs a moment to lock, and stale pre-retune samples
      sit in buffers);
   b. read exactly 15 s of samples = 36 million complex numbers. Each
      sample is 2 bytes (I and Q — think of them as the two coordinates
      of the radio wave's instantaneous state);
   c. chop them into 8,789 blocks of 4096, multiply each block by a
      window (tapers the edges so channels don't bleed), FFT each block
      (turns 4096 time samples into 4096 frequency channels), square
      (voltage→power), and **average all 8,789 spectra into one**. That
      average is the radiometer equation doing its work: the stored ON
      spectrum already has √8789 ≈ 94× less noise than a single FFT;
   d. same at 1423.2 MHz (OFF);
   e. append (timestamp, ON, OFF) to a buffer; every 20 cycles write the
      buffer to `chunk_NNNNNN.npz` — written to a temp name then renamed,
      so a power cut can never leave a half-written file.
4. **The log line** is your instrument panel:
   `cycle 42 14:23:05Z P_on/P_off=0.9710 P=+2.16 dB (31.0s)`
   — `P_on/P_off` near 1.000 (broadband power should be nearly equal at
   the two tunings; persistent deviation = something frequency-dependent,
   like the Sun's continuum, which you saw); `P` = total power relative
   to the session start (Sun arriving: smooth +2–3 dB; RFI: spiky;
   thermal drift: slow ±0.5 dB — all normal to see).
5. **Stopping** (Ctrl+C or systemd stop) flushes the buffer, finalizes
   `meta.json` (which records *everything*: frequencies, gain, pointing,
   site — so the session is self-describing forever).

What it does NOT do: no calibration, no science. It records. Everything
else happens later and can be redone with different settings any number
of times, which is why the raw chunks are sacred and everything else is
disposable.

## 2.2 `hi-process` — the refinery

Stages, in order, with what can go wrong at each:

1. **Load** all chunks, time-sorted.
2. **Quotient calibration** (§1.5 way 2). The OFF in the denominator is
   smoothed first — dividing by a *noisy* OFF would inject its noise into
   every channel; smoothing keeps its shape but not its noise. Output:
   spectra in Kelvin (via the assumed T_sys — a placeholder scale until
   hi-compare measures the real one).
3. **RFI flagging**, three passes, all in *channel* space because
   terrestrial interference lives at fixed frequency: (a) whole cycles
   with wild broadband power (arcing thermostat, garage door); (b)
   channels that are noisy across the whole session (a persistent
   carrier); (c) individual cycle×channel spikes. The velocity window
   where galactic HI can live is *protected* — otherwise the flagger
   would flag the Galaxy, which is, statistically speaking, a very
   suspicious signal. Watch the "flagged X%" line: 0–2% is normal; 10%+
   means a bad RFI environment — look at the waterfall to see what.
4. **Velocity regridding.** Each cycle gets its LSR correction (Earth's
   spin ±0.3 km/s, orbit ±30 km/s seasonal, solar motion 20 km/s — all
   computed by astropy from your site and the pointing) and is
   interpolated onto one common velocity grid. Side effect: any RFI that
   survived flagging smears out (it's fixed in frequency, so it lands in
   different velocity bins each cycle) while the sky stacks coherently.
5. **Baseline removal.** After calibration, slow curvature remains (the
   bandpass isn't *perfectly* identical 3.1 MHz apart). A low-order
   polynomial is fitted to each spectrum with asymmetric clipping —
   points well above the fit are excluded (that's emission!), points
   below only mildly — and subtracted. **This is the step that eats
   broad, weak emission if you're unlucky** (you saw it eat the
   inter-arm plateau in the PICTOR data); it's the known cost of
   quotient-style calibration, reduced by upgrade #2.
6. **Maps**: full time×velocity waterfall (`calibrated.npz`) and cycles
   averaged into 1° RA bins (`ra_map.npz`) — remember, in a drift scan
   *time is RA*: the sky slides through the beam at 15°·cos(dec)/hour.
7. **Extraction**: per RA bin — noise, peak, moments, N_HI column
   density, Gaussian components; globally — kinematic distances and the
   face-on map where geometry allows.

**How to read the plots** (`products/`):
* `waterfall_time.png` — time vs velocity. Healthy: smooth bright ridges
  drifting slowly. Vertical stripes = leftover RFI channels. Horizontal
  stripes = bad cycles. Everything black = check gain/tee.
* `avg_spectrum.png` — session average. Healthy: flat noise floor at 0 K
  with emission between roughly −100 and +50. Waves in the floor =
  baseline order too low.
* `waterfall_ra.png` — the science: your slice of the Milky Way.
* `coverage.png` — where that slice sits on the sky.

## 2.3 The rest, briefly but honestly

* **`hi-check`** = stethoscope. Run whenever anything about the hardware
  changed. Its three verdicts map to: ADC (is the digitizer in its happy
  range), bias tee (is the LNA alive), quotient noise vs radiometer
  limit (is the *whole chain* behaving like ideal physics — the single
  most powerful number; 1.0–2.0× is good).
* **`hi-sunscan`** = compass. The Sun is a ~70 K blowtorch: peak on it
  live (the P column), or let it drift through and fit the bump. Gives
  transit time → pointing, bump width → your real beam FWHM.
* **`hi-compare`** = the precision instrument. Simulates what HI4PI says
  your telescope should have seen for a grid of pointing errors,
  correlates each against what you actually saw, and reports the true
  pointing (typically to <1°, far finer than the beam) plus the true
  T_sys. Needs a processed session of ≥ several hours. The search is
  deliberately wide (±25°) — it assumes your pointing is worse than you
  think, which for a first-build wobbly mount is the correct assumption.
  If the correlation map's peak hugs the search edge, the true error is
  bigger still: rerun with a larger `compare.search_deg`.
* **`hi-stack`** = depth. Nights at the same dish setting average
  together; noise drops √N. `--auto` groups by declination so you don't
  have to remember which nights go together.
* **`hi-skymap`** = the progress bar.

**A decision table:**

| situation | command |
|---|---|
| touched any hardware | `hi-check` |
| moved the dish | sunscan (peak or transit), then update config |
| new night recorded | `hi-process --all` then `hi-stack --auto` |
| first processed night at a new pointing | `hi-compare` → put fitted az/el + T_sys in config → `hi-process --all --force` |
| curious how it's going mid-night | `hi-process <session> --last-h 2` |
| want to show someone | `hi-skymap`, `waterfall_ra.png`, `faceon_map.png` |

---

# Part III — Inspecting the physical build

The guiding principle: **losses before the SAWbird are catastrophic;
everything after it barely matters.** The SAWbird's 40 dB gain means the
rest of the chain contributes noise divided by 10,000. But anything lossy
*in front* of it attenuates the sky AND adds its own thermal noise. One
dB of loss before the LNA ≈ +60–70 K of T_sys — it can literally double
your noise. So spend your care budget in this order:

## 3.1 The probe-to-SAWbird path (criticality: 10/10)
The wire from the can's probe to the SAWbird input should be as short and
as direct as humanly possible — ideally the SAWbird bolted right at the
can with one connector. Every adapter ≈ 0.1–0.3 dB; every meter of cheap
coax ≈ 0.3–0.5 dB at 1.4 GHz; a corroded or loose connector can be 1+ dB.
Check: connectors metal-shiny, seated straight, snug (firm finger-tight
plus a nudge; don't gorilla them). If you have a jumper cable here,
making it shorter is the single highest-value hardware improvement
available to you.

## 3.2 The cantenna geometry (criticality: 7/10)
The can is a circular waveguide; it only works if the wave fits:
* **Diameter D**: must exceed 12.4 cm (below that, 1420 MHz physically
  cannot propagate — the can is dead). Ideal 15–18 cm. Bigger than ~20 cm
  starts supporting unwanted modes.
* **Probe length**: a quarter of the free-space wavelength, ≈ **53 mm**,
  measured from the can wall's inner surface. ±3 mm is fine.
* **Probe position**: a quarter of the *guide* wavelength from the closed
  back wall. The guide wavelength depends on D:

  | D (cm) | λg (cm) | probe at (cm from back) |
  |---|---|---|
  | 15 | 37.3 | 9.3 |
  | 16 | 33.3 | 8.3 |
  | 17 | 31.0 | 7.8 |
  | 18 | 29.0 | 7.3 |

* **Can length**: at least ~¾ λg (≈ 22–28 cm) from back wall to opening.
* Dents/ovality: a few mm is noise; a crushed can is not.
* Probe soldering: solid mechanical joint to the connector center pin;
  a cold joint here is a silent 3 dB.

Getting these to ±10% costs maybe 1 dB versus perfect. Real, not fatal.

## 3.3 Feed placement at the dish (criticality: 6/10 — and measurable!)
The can's opening should sit where the old LNB's feed horn sat, opening
facing the *center of the dish* (for an offset dish that means the can
looks *up* at the dish at an angle — this is why the beam points ~20–25°
above where the dish face seems to aim, your "offset" question).
You do not need to compute anything: **use the Sun as your optimizer.**
With `hi-observe` running and the beam peaked on the Sun, slide the can
a few cm in/out along its axis and watch the `P` column — maximize it.
Then tilt the can a degree or two — maximize again. That empirically
finds focus better than any tape measure. Probe rotation, by the way,
does not matter for HI (the line is unpolarized) — don't waste time on it.

## 3.4 The mount wobble (criticality: 4/10 — less than you fear)
Wobble smears your beam: a ±2° wander on a 17° beam is a ~3% effect —
invisible. What wobble actually costs you is *repeatability*: if the dish
sags overnight or shifts in wind, tonight's pointing isn't yesterday's.
Mitigations, all cheap: guy the mount with string/wire after aiming; mark
positions (pencil lines, zip ties) so settings are reproducible; and let
software absorb the rest — `hi-compare` measures each session's true
pointing after the fact, so even a slow sag becomes a known number rather
than an error. A wobbly mount with per-session hi-compare beats a rigid
mount with assumed pointing.

## 3.5 Ground pickup and surroundings (criticality: 5/10)
The feed's backside and the dish's spillover "see" the ground, which
glows at ~300 K — this is usually the biggest avoidable T_sys
contribution in amateur builds. Rules of thumb: keep elevation above
~30° where possible; avoid pointing over buildings/metal roofs; a simple
aluminum-foil skirt around the can's opening (a crude choke) can shave
tens of K. Trees are surprisingly benign; warm brick walls are not.

## 3.6 Weather (criticality: 2/10 for physics, yours to judge for safety)
Rain on the *sky path* is irrelevant at 21 cm — this band laughs at
clouds. Water *on the feed* is not: a wet probe or a puddle in the can
detunes and attenuates. A plastic bag or thin tarp over the can
(polyethylene is transparent at 21 cm) plus a drain hole in the can's
low point handles dew and drizzle. Electronics in a food container with
the coax exiting downward (drip loop). Your "don't observe in rain"
policy is completely fine — HI accumulates whenever you do run, and the
stacking doesn't care about gaps.

## 3.7 What "bad" looks like in the data (your build-quality meter)
You don't need to *guess* build quality — measure it:
* `hi-check` quotient noise vs radiometer limit: ≤2× = chain is fine
  regardless of how it looks.
* Sunscan bump height: remeasure after any hardware change; if a change
  raises the Sun bump, it lowered your losses. The Sun is your free
  signal generator.
* `hi-compare` gain → implied T_sys: THE scoreboard. 60–90 K excellent,
  90–150 K typical first build, 150–250 K = go hunting front-end losses
  (§3.1), >300 K = something is broken, not merely imperfect.

---

# Part IV — The three upgrades, explained properly

## 4.1 Continuum products (easy, immediate)
HI is a *line* — one narrow frequency. But the sky also glows across
*all* frequencies (synchrotron electrons, HII regions, the Sun, Cas A):
the *continuum*. Your chunks already contain it — the broadband power of
every cycle (the sunscan uses it; nothing else does yet). The upgrade:
hi-process additionally outputs power-vs-RA per session, giving you a
continuum drift profile alongside every HI map. What you'd see: the
galactic plane as a broad bump, the Sun if it's up, and plausibly Cas A
(the sky's second-brightest radio source, dec +59° — beautifully placed
for Lithuania) as a repeatable ~0.5–1 K transit spike. Why care beyond
the extra science: a repeatable point-source transit is a *daily free
calibrator* — pointing and gain checked every day without touching
anything. Cost: modest additions to process.py; the data already exists.

## 4.2 In-band frequency switching + folding (the sensitivity upgrade)
Today: ON has the line, OFF (3.1 MHz away) has nothing — 50% of time
spent measuring pure calibration. The upgrade: make the two tunings only
~0.7 MHz apart, so **the line is inside the band at both** — just landing
in different channels (shifted by 0.7 MHz ≈ 148 km/s). Form the quotient
as usual: the bandpass still cancels (it's still glued to the LO), but
now the result contains the line *twice*: once positive (from ON) and
once negative, shifted (from OFF). Then **fold**: shift the negative copy
back by +0.7 MHz, flip its sign, and average it with the positive copy.
Both halves of every cycle now contribute line signal → noise improves by
√2 → **+41% sensitivity, or equivalently every night counts as two**.
Bonus: baselines get *better*, because the two tunings are 0.7 MHz apart
instead of 3.1 — their bandpasses are more nearly identical, so less
residual curvature for the polynomial to absorb (the thing that ate the
PICTOR plateau). Subtleties to handle (why it's a weekend, not an
evening): the shift must exceed the line's full velocity extent or the
line partially cancels itself (0.7 MHz ≈ 148 km/s is safe for
|b| > few°, marginal dead on the plane — make it configurable); the DC
spike lands at a different velocity in each tuning (blank both); and the
folding must happen after LSR regridding. Plan: implement behind a
config flag (`switching: inband` vs `classic`), verify with
`hi-simulate` (inject known sky, check amplitude and √2 noise), then one
real A/B night. Old sessions remain processable — the mode is recorded
in meta.json.

## 4.3 CI smoke test (insurance)
A GitHub Actions workflow that, on every push, runs the simulator (short
session, known injected sky), processes it, and asserts: products exist,
the recovered line amplitude matches the injected one to tolerance, and
the noise matches the radiometer prediction. Why bother when you can test
by hand: because the failure mode of *calibration* code is silent
wrongness — a sign flip or an off-by-one channel doesn't crash, it just
quietly makes every map subtly wrong. The end-to-end sim catches exactly
this class (it already caught real bugs during development); CI makes it
impossible to forget to run. Cost: one YAML file and a small assertion
script; free on public repos.
