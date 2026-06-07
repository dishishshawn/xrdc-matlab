# Validation results — Tushar

Sample **S25 — SRO on STO(100)**, 700 °C / 100 mT / 10500 sh / 5 Hz (2026-05-06).
App run via `xrdc.io.readScan` → `findPeaks` / `fitPeak` / `analyzeFringes`.
All three figures reproduced by `runValidation.m` (2026-06-06); θ-2θ first run 2026-06-05.

| Date       | Input file              | Scan type | Key metric (app)                          | Tushar          | Match? | Notes |
| ---------- | ----------------------- | --------- | ----------------------------------------- | --------------- | ------ | ----- |
| 2026-06-05 | `…_2theta_omega_…txt`   | θ-2θ      | STO 100/200/300 = 22.766 / 46.493 / 72.595° | graph (same)    | ✓      | Matches theory (22.754/46.472/72.567°) to <0.03°; film shoulders + noise floor match his figure. `output/S25_2theta_omega_xrdc.png`. |
| 2026-06-05 | `…_sub RC_…txt`         | ω (RC)    | FWHM 0.0304° **(gauss)** / 0.0253° (lorentz) | 0.0308°         | ✓*     | Gaussian fit matches to −1.3%; app's default Lorentzian reads −18% low. |
| 2026-06-05 | `…_film RC_…txt`        | ω (RC)    | FWHM 0.0911° **(gauss)** / 0.0776° (lorentz) | 0.0920°         | ✓*     | Gaussian fit matches to −1.0%; Lorentzian −15.7% low. Combined RC figure: `output/S25_subfilmRC_xrdc.png`. |
| 2026-06-06 | `…_XRR_…_11.txt`        | XRR       | d = **40.3 ± 0.4 nm** (21 fringes, quad Kiessig) | graph (fringes) | ✓      | XRR `.txt` now supplied. Fringe count, spacing, TER edge (~0.61°) and decay envelope reproduce his figure; quad/linear/FFT thicknesses agree (40.3/41.6/40.7 nm). `output/S25_XRR_xrdc.png`. |

`✓* = matches once fit shape is set to Gaussian.`

## Headline
- **θ-2θ and core I/O: validated.** Peak positions agree with Tushar and with first-principles
  STO spacings to better than 0.03°. `output/S25_2theta_omega_xrdc.png` reproduces his plot
  (peaks, film shoulders, log-Y noise floor).
- **Rocking-curve FWHM: matches to ~1% — but only with `Shape="gauss"`.** Tushar evidently fit
  Gaussians; the app defaults to **Lorentzian**, which sits ~16–18% low on both curves. This is a
  fit-model choice, not a data or algorithm error (centres and peak detection agree). The combined
  normalised figure `output/S25_subfilmRC_xrdc.png` overlays both curves the way his does
  (sub black / film red, Δω axis, FWHM in legend).
- **XRR: validated.** With the XRR `.txt` now in hand, `analyzeFringes` finds 21 Kiessig fringes
  and returns d = 40.3 ± 0.4 nm, with linear-sinθ (41.6 nm) and FFT (40.7 nm) cross-checks within
  ~3%. `output/S25_XRR_xrdc.png` reproduces his direct-beam drop, TER plateau peak, fringe
  spacing, and noise floor.

## Fit-shape sweep (FWHM, degrees)
| Curve   | Tushar | lorentz (default) | gauss   | pseudoVoigt | half-max |
| ------- | ------ | ----------------- | ------- | ----------- | -------- |
| sub RC  | 0.0308 | 0.0253 (−18.0%)   | 0.0304 (−1.2%) | 0.0274 (−10.9%) | 0.0200 (−35%) |
| film RC | 0.0920 | 0.0776 (−15.7%)   | 0.0911 (−1.0%) | 0.0892 (−3.1%)  | 0.0850 (−7.6%) |

## Reproduce
From the repo root, in MATLAB: `run('validation/tushar/runValidation.m')`. It regenerates the
XRR and combined-RC PNGs in `output/` and prints the FWHM / thickness summary.

## Follow-ups
1. Confirm with Tushar that his FWHM came from a **Gaussian** fit. If so, document it and
   consider whether RC analysis should default to (or offer) Gaussian.
2. ~~Ask Tushar for the XRR `.txt`.~~ **Done** — XRR validated, d = 40.3 nm. Optionally ask Tushar
   for *his* fitted thickness to put a number against ours (his figure has no thickness annotation).
