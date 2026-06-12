# Design: Auto min-prominence for θ-2θ peak detection

**Date:** 2026-06-12
**Status:** Approved (Shawn, via design review in session)
**Feature:** `xrdc.peaks.findPeaks` auto-selects a peak-prominence criterion that finds
the physically correct peak set on θ-2θ scans without manual threshold tuning; the GUI
θ-2θ panel defaults to it.

## Problem

θ-2θ intensity spans ~6 decades: substrate peaks ~10⁶ counts, film (00l) peaks
~10²–10³, background ~10–100. No single linear prominence threshold works — high
enough to reject noise it misses film peaks; low enough for film peaks it admits
noise and over-segments.

Both failure modes are live today:

- The GUI default (`MinProminence = 5% of max counts`) misses film peaks entirely on
  Tushar's PTO/STO scans; he manually lowers "Min prom (%)" per scan until the right
  peaks appear.
- The library default (1.5% of data range) over-segments: DATA_SWEEP #16/#7 report
  66–69 "peaks" on S10/S06, with the Kα-split substrate peak counted several times
  (DATA_SWEEP finding #3).

Ground truth for "correct" on a PTO/STO θ-2θ scan: substrate (00l) series + film
(00l) series only. No Kα₁/Kα₂ duplicates, no Laue-fringe wiggles, no noise spikes.

## Approach (chosen: log-domain detection)

Humans judge peaks on the log plot: a peak is real when it stands out by a fraction
of a decade above its surroundings, regardless of absolute counts. The auto path
detects on `log10(counts)` with a prominence threshold in **decades** — scale-invariant
across the full dynamic range, so one default covers substrate and film peaks.

Rejected alternatives:

- **Adaptive linear threshold (k × MAD noise)** — minimal change but does not fix
  over-segmentation: fringe wiggles and shoulders on strong peaks have prominence ≫
  noise; the S06/S10 blowups survive.
- **Plateau sweep (auto-tune until peak count stabilizes)** — on STO/PTO scans the
  widest stable plateau is "substrate only" (film peaks drop out around 10³ counts
  while substrate peaks survive to 10⁶), so it converges on the wrong answer for
  exactly the target data. Slower and unpredictable as a default.

## Components

### 1. `xrdc.peaks.findPeaks` — auto path (modified)

Triggered exactly as today: `MinProminence` omitted / NaN. Replaces the
`max(1, 0.015·range)` rule (no opt-in flag; the old rule is removed). Explicit numeric
`MinProminence` → the existing linear path, behaviour unchanged.

Auto pipeline (all within the `TwoThetaRange` window, as today):

1. **Log transform.** `yLog = log10(max(y, 1))`. log10 is monotone, so local-maxima
   indices on `yLog` coincide exactly with those on `y` — detect on log, report from
   linear, no index mismatch.
2. **Log-domain detection.** `findpeaks(yLog)` with
   `MinPeakProminence = 0.3` (decades — peak stands ~2× above its neighbouring
   troughs) and the `MinPeakDistance` derived from `MinSeparation` as today.
3. **Noise-floor guard (Poisson significance).** Reject candidates with
   `y(pk) − bg(pk) < 5·sqrt(max(bg(pk), 1))` where `bg = movmedian(y, w)` and `w`
   spans ~1° of 2θ (converted to samples via the median step). This kills the large
   fake log-prominence of 1-vs-3-count quantisation jitter at the baseline — the
   known failure mode of pure log detection.
4. **Linear metrics.** A second `findpeaks` pass on linear `y` with near-zero
   prominence returns every local maximum with its linear FWHM / prominence / height;
   intersect by index with the accepted log-domain set. `MinHeight`, `MinWidth`,
   `MaxWidth` are enforced on these linear values, same semantics as the manual path.

Returned struct array is unchanged in shape and semantics (counts, prominence, fwhm,
leftHalf, rightHalf, index — all linear-domain). Downstream consumers
(`fitPeak`, `adjustPeaks`, `identifyMaterial`) see no interface change.

**Constants.** 0.3 decades, 5σ, and the ~1° background window are internal constants,
documented in the function help and SCIENTIFIC_ASSUMPTIONS — not API options (add
knobs only if real use demands them). Initial values are educated guesses; they are
**calibrated against the seven Tushar θ-2θ scans in `data/` during implementation**
(explicit plan task). The acceptance criteria below are the commitment; the constants
may move to meet them.

**Fallback path.** When the Signal Processing Toolbox is absent,
`xrdc.peaks.findpeaks_fallback` is used for both passes, as today. The fallback must
support the subset of options the auto path uses (verify during implementation; extend
it if a used option is missing).

### 2. GUI — `xrdcApp.m` θ-2θ panel (modified)

- "Min prom (%)" edit field default becomes the string `auto`.
  - Value `auto` (case-insensitive) or empty → call `findPeaks` with **no**
    `MinProminence` and `MinSeparation = 0.2°`. The separation bump merges the
    Kα₁/Kα₂ substrate split (~0.12° at the STO 002) per DATA_SWEEP finding #3.
  - Numeric value → exactly today's manual behaviour:
    `MinProminence = max(counts)·pct/100`, library-default separation. Tushar's
    existing habit keeps working.
  - Anything else → treated as `auto` (forgiving parse), no error dialog.
- `runThetaTwoTheta` and `onIdentifyMaterial` currently duplicate the peak-detection
  call; extract one shared helper (e.g. `thetaPeaks(st)`) so the plot and Identify
  Material always operate on the same peak set.
- The "No peaks detected - lower Min prom" alert text updates to mention `auto`.
- Library `MinSeparation` default (0.05°) is untouched — superlattice and validation
  scripts pass explicit values and must not shift.

### 3. Docs

- `docs/USER_GUIDE.md`: θ-2θ parameter table ("Min prom (%)": default `auto`, what it
  does, how to override); workflow text that says "lower the threshold" updated.
- `docs/SCIENTIFIC_ASSUMPTIONS.md`: new entry — log-domain prominence algorithm,
  constants, Poisson guard, and the divergence note (the Delphi original has no auto
  mode; slope-threshold detectors only).
- `docs/FEATURES.md`: status line.
- DATA_SWEEP finding #3 gets a "fixed by auto prominence" follow-up note.

## Error handling

- Negative explicit `MinProminence` → existing `xrdc:peaks:badProminence` error
  (unchanged).
- All-zero / flat scans → auto path returns the empty peak array (guard step rejects
  everything), no error.
- Scans with < 3 points, missing fields, size mismatch → existing errors unchanged.

## Testing

`tests/testPeaks.m` additions (functiontests, per repo conventions):

- **Synthetic known-answer:** 10⁶-count substrate Gaussian + 10³-count film Gaussian
  on ~30-count Poisson-like background (fixed rng) → auto finds exactly 2 peaks at
  the right positions.
- **Noise guard:** flat quantised background (counts ∈ {0..3}) → 0 peaks.
- **Fringe suppression:** film peak with superimposed small thickness-fringe ripple →
  fringes below 0.3 decades are not reported.
- **Backward compatibility:** explicit `MinProminence` on a synthetic scan returns
  results identical to the pre-change behaviour (pin current output).
- **Option interplay:** `MinHeight` / `TwoThetaRange` / `MinSeparation` respected in
  auto mode.
- **Real data (gated `assumeTrue(isfile(...))`, never required for CI green):**
  - S04, S05, S11 θ-2θ → substrate STO 001/002/003 within 0.05° of 22.75/46.47/72.57°
    plus film peaks; exact expected number of peaks per file pinned after the
    calibration task.
  - S10 (the 66-peak blowup) → no two reported peaks within 0.15° of each other;
    fewer than 20 peaks total.
  - S06 (PTO/LAO, the 69-peak case) → same over-segmentation bound.

GUI helper logic (`auto` vs numeric parse) covered by a unit test if the helper is
factored to be callable headless; otherwise verified in the integration review.

## Out of scope (v1)

- Auto-prominence for XRR fringe analysis (`analyzeFringes` has its own log-decade
  prominence already) and φ scans (`findPhiPeaks` has a noise-σ criterion).
- Exposing decade-threshold / guard constants as API or GUI options.
- Kα₂ stripping / profile deconvolution (separation merge is sufficient for ID).
- Regenerating the DATA_SWEEP table (ad-hoc script, re-run any time).
