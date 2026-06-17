# Feature checklist — xrdc-matlab

Status of every feature, derived from the actual `+xrdc/` source (not the roadmap).
Legend: ✅ done · ⏳ in progress / gated · ⬜ wanted, not started · ❌ stubbed (raises error).

_Last reconciled: 2026-06-09._

## File I/O (`+xrdc/+io`)
- ✅ PANalytical `.xrdml`, single scan — `readXrdml`
- ✅ PANalytical `.xrdml`, multi-scan / area (RSM) — `readXrdmlArea`
- ✅ Rigaku SmartLab `.txt`, headered **and** headerless — `readRigakuTxt` / `readTextScan`
- ✅ Rigaku SmartLab Studio II `.hgx` (GlobalFit, HDF5) — `readRigakuHgx` (curve verified bit-identical to the `.txt` twin)
- ✅ Philips `.x00` — `readPhilipsX00`
- ✅ Plain two-column text — `readTextScan`
- ✅ Auto-detect format from first 512 bytes — `readScan`
- ✅ German-locale comma-decimal tolerance on input
- ✅ Common scan struct shape — `emptyScan`
- ❌ Rigaku binary `.raw` — `readRigakuRaw` raises `xrdc:io:notImplemented` (awaiting sample files)
- ❌ Rigaku ASCII `.ras` — `readRigakuRas` raises `xrdc:io:notImplemented` (awaiting sample files)
- ✅ Rigaku `.hgx` (HDF5) — reader added and tested; layer-model parameters under
      `/current/parameters` not yet parsed (curve + simulated fit are)

## Scan types & analyses
- ✅ θ-2θ (`twoThetaOmega`) — peaks + substrate overlay **[validated vs Tushar S25]**
- ✅ Material identification from (00l) peaks — `xrdc.lattice.identifyMaterial`: harmonic-series grouping, Kβ/W-Lα ghost filter, strain-aware candidate matching (pseudomorphic c prediction), ranked candidate sets, PZT composition estimate (strain-caveated). Database: `+xrdc/+data/materials.json` (STO, SRO, PTO, PZT). **[validated vs S25 + S31]**
- ✅ Rocking curve (`omega`) — peak fit → FWHM **[validated vs Tushar S25, Gaussian]**
- ✅ φ scan (`phi`) — noise-floor pole detection + 360° wrap handling + n-fold
      symmetry report (`findPhiPeaks`); robust on weak films **[validated on PtO₂/TiO₂(101)]**
- ✅ XRR — critical edge, Kiessig fringes, film thickness (`+xrr`)
- ✅ XRR slab-model fit — Parratt + Névot–Croce, fits thickness/density/roughness with Jacobian errors (`fitReflectivity`, `reflectivityModel`, `opticalConstants`)
- ✅ Reciprocal-space map (`area`) — load, transform, contour (`+rsm`)
- ✅ **RSM strain & composition** (`xrdc.rsm.analyzeStrainRSM`) — DONE. Auto-detects
  substrate+film peaks in an asymmetric RSM, recovers a∥/a⊥, biaxial strain,
  relaxed a0, relaxation, and PZT composition for a declared film. Geometry
  validated on real TiO₂ substrate; strain/composition synthetic-only this
  iteration (see SCIENTIFIC_ASSUMPTIONS). Candidate ranking / identifyMaterial
  integration is a deliberate follow-on.
- ✅ Scan-type auto-detection from file content + name

## Peak detection & fitting (`+xrdc/+peaks`)
- ✅ Peak detection — `findPeaks` (Signal Processing, with `findpeaks_fallback`)
- ✅ Auto prominence — log-domain (0.3 decades) + unit-aware Poisson guard; default when `MinProminence` omitted
- ✅ Delphi-compatible legacy mode — `findPeaksLegacy`
- ✅ Profile fits — `fitPeak`, shapes: **lorentz** (default) · **gauss** · **pseudoVoigt**
- ✅ Jacobian-based parameter standard errors (Optimization Toolbox path)
- ✅ Multistart robust fitting (Global Optimization path)
- ✅ Manual peak adjustment — `adjustPeaks`
- ✅ GUI rocking-curve panel now reports FWHM under the **other** fit shapes as a cross-check
      (surfaces the ~16% Gaussian-vs-Lorentzian gap instead of hiding it behind the default)
- ⬜ Decision still open: should the RC *default* be Gaussian? (Tushar reports Gaussian; app
      default is still Lorentzian — left as a deliberate decision, not changed silently)

## Lattice / crystallography (`+xrdc/+lattice`)
- ✅ d-spacing from (hkl) per crystal system — `dSpacingFromHKL`
- ✅ 2θ ↔ d conversions — `dToTwoTheta`, `twoThetaToD`
- ✅ Energy ↔ wavelength — `energyToLambda`, `lambdaToEnergy`
- ✅ Nelson–Riley lattice-parameter refinement — `nelsonRiley`
- ✅ Bragg-pattern simulation — `simulatePattern`
- ✅ Kiessig thickness from fringe spacing — `thicknessFromFringes`
- ✅ Superlattice period Λ from satellite spacing — `superlatticePeriod`
      (per-layer A/B thickness needs Parratt/optical modelling — out of scope)

## Signal processing (`+xrdc/+signal`)
- ✅ Smoothing — `smoothCounts` (csaps spline / sgolay / `sgolay_fallback`)
- ✅ Background subtraction — `subtractBackground`
- ✅ Savitzky–Golay derivatives — `derivatives`
- ✅ XRR log-envelope — `logEnvelope`

## Plotting (`+xrdc/+plot`)
- ✅ Single-scan publication plot — `plotScan` (log-Y, peak marking)
- ✅ RSM filled contour w/ log-decade colorbar — `plotRsm`
- ✅ Stacked multi-scan — `plotStack`
- ✅ Schwaigert et al. Fig 2 styling — `publicationStyle`
- ✅ 600 dpi export

## GUI (`xrdcApp.m`)
- ✅ Load → auto-detect → analyze → live re-plot → 600 dpi export
- ✅ Zero-MATLAB-knowledge workflow for lab members
- ✅ Animated launch splash (`xrdc.ui.showSplash`, `resources/splash.html`)
- ✅ "Customize Plot…" — publication overrides on top of defaults
      (`xrdc.plot.applyStyle`): title/labels, axis limits, log/linear, font
      size, line colour/width, markers, grid, and sized export (in)
- ✅ "Settings…" — lab presets persisted as JSON in `prefdir`
      (`+xrdc/+config`: `defaults`/`load`/`save`/`configPath`). Overrides four
      baked-in Paik-lab defaults without code edits: X-ray wavelength (stamped
      onto formats that drop it — Rigaku `.txt`/`.hgx`, plain text; `.xrdml`/
      `.x00` keep their real λ), default substrate (material ID), rocking-curve
      default fit shape, and the filename→scan-type rules (toggle + fallback
      scan type). Unknown/corrupt settings degrade to defaults.

## Toolbox fallbacks (run without licenses)
- ✅ Signal Processing fallback (`findpeaks_fallback`, `sgolay_fallback`)
- ✅ Optimization fallback (brute-force fit when no `lsqcurvefit`)
- ✅ Curve Fitting optional (`csaps` when present)

## Testing
- ✅ Suites for io, lattice, peaks, plot, rsm, signal, xrr (`tests/test*.m`)
- ✅ Full `runtests` green on the lab MATLAB (R2026a) — 105 pass, 0 fail, 5 skipped
      (real-data/toolbox-gated). Handoff gate met.

## Distribution / packaging
- ✅ MATLAB-project usage (clone + `addpath`)
- ✅ Standalone `XRDC.exe` — builds via `build/buildStandalone.m` (Compiler R2026a);
      single-file exe produced. Defender exclusion on the build folder needed for the
      embed step, or use `buildStandalone(SingleFile=false)`. Runtime-verify on a
      clean (no-MATLAB) machine still pending.

## Validation against external references
- ✅ θ-2θ — Tushar S25 SRO/STO(100): peak positions match to <0.03°
- ✅ Rocking-curve FWHM — Tushar S25: matches to ~1% with Gaussian fit
- ✅ XRR / Kiessig thickness — S25 XRR (found in `data/`) → 40.30 ± 0.44 nm (21 fringes);
      full `data/` sweep (37 files, 0 crashes) in `validation/DATA_SWEEP.md`
- ✅ XRR — our ~40 nm (fringe spacing) vs Tushar's GenX fit ~36 nm reconciled: not a bug.
      Thickness is flat vs low-angle cutoff (refraction ruled out); gap is the expected
      fringe-spacing-vs-Parratt difference. Only GenX model params from Tushar would close it.
- ⬜ RSM strain/composition — geometry (hkl inversion) validated on TiO₂ substrate of real PtO₂/TiO₂ 112 map (a=4.593, c=2.959 Å); biaxial strain/composition validated synthetic-only this iteration

## Planned / wanted (from README, not started)
- ⬜ Statistics & ML: bootstrap CIs on fit parameters, outlier rejection on peak picks
- ⬜ Image Processing: 2D detector images (GIWAXS / RSM area images)
- ⬜ MATLAB Compiler standalone distribution (see Packaging above)
- ✅ Lab presets / settings menu — DONE. Several defaults were baked in for the Paik
      lab's Rigaku and would be wrong at another lab. The `xrdcApp` "Settings…" menu now
      overrides these without code edits, persisted as JSON in `prefdir` (`+xrdc/+config`):
    - **X-ray wavelength** — readers still hardcode `lambda = 1.5406` (Cu Kα₁) where the
      `.txt`/`.hgx` export drops it, but the app re-stamps the configured λ onto those
      scans on load (formats that carry a real λ — `.xrdml`/`.x00` — are left untouched).
    - **Default substrate** — the app's material-ID substrate dropdown defaults to the
      configured substrate (was hardcoded "SrTiO3"). `analyzeStrainRSM`'s own default is
      unchanged (no GUI entry point yet).
    - **Filename → scan-type rules** — the `TR_` prefix + `RC`/`XRR`/`2theta` cues can be
      switched off (`useFilenameRules`) for labs with a different naming convention, with
      a configurable fallback scan type when nothing else disambiguates.
    - **Rocking-curve default fit shape** — the app's RC shape dropdown default is now
      configurable (gauss / lorentz / pseudoVoigt); resolves the open Gaussian-vs-Lorentzian
      default by making it a lab setting.

---
_Add new "wanted" items below as they come up so this stays the single source of truth._
