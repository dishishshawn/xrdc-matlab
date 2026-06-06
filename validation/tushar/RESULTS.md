# Validation results — Tushar

Sample **S25 — SRO on STO(100)**, 700 °C / 100 mT / 10500 sh / 5 Hz (2026-05-06).
App run via `xrdc.io.readScan` → `findPeaks` / `fitPeak` on 2026-06-05.

| Date       | Input file              | Scan type | Key metric (app)                          | Tushar          | Match? | Notes |
| ---------- | ----------------------- | --------- | ----------------------------------------- | --------------- | ------ | ----- |
| 2026-06-05 | `…_2theta_omega_…txt`   | θ-2θ      | STO 100/200/300 = 22.766 / 46.493 / 72.595° | graph (same)    | ✓      | Matches theory (22.754/46.472/72.567°) to <0.03°; film shoulders + noise floor match his figure. |
| 2026-06-05 | `…_sub RC_…txt`         | ω (RC)    | FWHM 0.0304° **(gauss)** / 0.0253° (lorentz) | 0.0308°         | ✓*     | Gaussian fit matches to −1.2%; app's default Lorentzian reads −18% low. |
| 2026-06-05 | `…_film RC_…txt`        | ω (RC)    | FWHM 0.0911° **(gauss)** / 0.0776° (lorentz) | 0.0920°         | ✓*     | Gaussian fit matches to −1.0%; Lorentzian −15.7% low. |
| —          | (no XRR `.txt` supplied) | XRR       | —                                         | graph present   | ✗ data | Reference XRR figure is in `reference/`, but no XRR scan in `input/` — can't reproduce yet. |

`✓* = matches once fit shape is set to Gaussian.`

## Headline
- **θ-2θ and core I/O: validated.** Peak positions agree with Tushar and with first-principles
  STO spacings to better than 0.03°. The exported `output/S25_2theta_omega_xrdc.png` reproduces
  his plot (peaks, film shoulders, log-Y noise floor).
- **Rocking-curve FWHM: matches to ~1% — but only with `Shape="gauss"`.** Tushar evidently fit
  Gaussians; the app defaults to **Lorentzian**, which sits ~16–18% low on both curves. This is a
  fit-model choice, not a data or algorithm error (centres and peak detection agree).

## Fit-shape sweep (FWHM, degrees)
| Curve   | Tushar | lorentz (default) | gauss   | pseudoVoigt | half-max |
| ------- | ------ | ----------------- | ------- | ----------- | -------- |
| sub RC  | 0.0308 | 0.0253 (−18.0%)   | 0.0304 (−1.2%) | 0.0274 (−10.9%) | 0.0200 (−35%) |
| film RC | 0.0920 | 0.0776 (−15.7%)   | 0.0911 (−1.0%) | 0.0892 (−3.1%)  | 0.0850 (−7.6%) |

## Follow-ups
1. Confirm with Tushar that his FWHM came from a **Gaussian** fit. If so, document it and
   consider whether RC analysis should default to (or offer) Gaussian.
2. Ask Tushar for the **XRR `.txt`** matching the XRR reference figure to validate the Kiessig
   thickness path (the one analysis his graph set covers that we can't yet check).
