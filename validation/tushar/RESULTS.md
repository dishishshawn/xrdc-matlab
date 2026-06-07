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

## Tushar confirmations (2026-06-07, email)
- **RC fit shape = Gaussian.** Tushar confirmed he used a Gaussian fit for the FWHM values.
  This matches our finding exactly. The GUI rocking-curve panel now **defaults to Gaussian**
  (`xrdcApp.m`, `case 'omega'`); Lorentzian/pseudoVoigt still selectable, and all shapes are
  shown as a cross-check.
- **XRR thickness (GenX) ≈ 36 nm.** Our Kiessig fringe-spacing estimate is 40.3 nm (quad), with
  linear 41.6 and FFT 40.7 — internally consistent at ~40–41 nm, i.e. **~11% above** his GenX
  value. GenX fits the full reflectivity curve (critical edge + roughness + density), so some
  divergence from a pure fringe-spacing estimate is expected; worth reconciling but not a bug.

## XRR thickness reconciliation (2026-06-07) — not a bug
Investigated the ~4 nm gap (our ~40 nm fringe-spacing vs Tushar's GenX ~36 nm).
**Conclusion: our estimate is sound; the gap is a method difference, not an error.**
- **Refraction/critical-edge ruled out.** Re-fitting with only high-angle fringes (where
  refraction is negligible) leaves the thickness flat: linear fit = 41.6 (all) → 41.1 (>1°) →
  40.6 (>1.5°) → 40.4 (>2°) → 40.8 (>2.5°). A critical-edge-handling bug would pull the
  high-angle value down toward 36; it doesn't. Evidence: `output/S25_XRR_reconcile.png`.
- **What's left** is the expected gap between a geometric fringe-spacing estimate and a full
  Parratt fit (GenX models roughness + density + possibly a 2-layer stack, which can legitimately
  yield a smaller *layer* thickness than the total-stack fringe period).
- **Minor quality note:** a couple of detected fringes are slightly off-cadence (~0.9° and ~4.3° —
  one likely split/missed fringe). The robust 21-point linear fit absorbs it; it adds scatter, not
  bias. Tightening fringe detection is optional and would not close the gap.
- **To fully close it** we'd need Tushar's GenX model (layer structure, roughness, density). Asking
  him is the only remaining step — code side is done.

## Follow-ups
1. ~~Confirm Gaussian RC fit.~~ **Done** — confirmed; RC default flipped to Gaussian.
2. ~~Ask Tushar for XRR thickness.~~ **Done** — GenX ≈ 36 nm vs our ~40 nm (see above).
3. ~~Reconcile the ~4 nm XRR gap.~~ **Done** — not a bug (see reconciliation above); only open
   piece is getting GenX model params from Tushar.
4. **Superlattice / heterostructure (Tushar request).** Implemented: `xrdc.lattice.superlatticePeriod`
   (period Λ from satellite spacing) + tests + `demoSuperlattice.m`, synthetic-validated to 0.01%.
   Still needs validation against a **real** superlattice θ-2θ scan — ask Tushar to send one.
   Per-layer thickness needs full Parratt modelling (GenX territory) — out of scope.
