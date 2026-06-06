# Feature checklist — xrdc-matlab

Status of every feature, derived from the actual `+xrdc/` source (not the roadmap).
Legend: ✅ done · ⏳ in progress / gated · ⬜ wanted, not started · ❌ stubbed (raises error).

_Last reconciled: 2026-06-05._

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
- ✅ Rocking curve (`omega`) — peak fit → FWHM **[validated vs Tushar S25, Gaussian]**
- ✅ φ scan (`phi`) — linear-Y, 4-fold symmetry detection
- ✅ XRR — critical edge, Kiessig fringes, film thickness (`+xrr`)
- ✅ Reciprocal-space map (`area`) — load, transform, contour (`+rsm`)
- ✅ Scan-type auto-detection from file content + name

## Peak detection & fitting (`+xrdc/+peaks`)
- ✅ Peak detection — `findPeaks` (Signal Processing, with `findpeaks_fallback`)
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

## Toolbox fallbacks (run without licenses)
- ✅ Signal Processing fallback (`findpeaks_fallback`, `sgolay_fallback`)
- ✅ Optimization fallback (brute-force fit when no `lsqcurvefit`)
- ✅ Curve Fitting optional (`csaps` when present)

## Testing
- ✅ Suites for io, lattice, peaks, plot, rsm, signal, xrr (`tests/test*.m`)
- ⏳ Full `runtests` green on the **lab** MATLAB — stated handoff gate (pending)

## Distribution / packaging
- ✅ MATLAB-project usage (clone + `addpath`)
- ⏳ Standalone `XRDC.exe` — `build/buildStandalone.m` + checklist scaffolded;
      compile gated last, pending MATLAB Compiler **install** (licensed, not installed)

## Validation against external references
- ✅ θ-2θ — Tushar S25 SRO/STO(100): peak positions match to <0.03°
- ✅ Rocking-curve FWHM — Tushar S25: matches to ~1% with Gaussian fit
- ✅ XRR / Kiessig thickness — S25 XRR (found in `data/`) → 40.30 ± 0.44 nm (21 fringes);
      full `data/` sweep (37 files, 0 crashes) in `validation/DATA_SWEEP.md`
- ⬜ XRR — still worth confirming the 40.3 nm against Tushar's reported value
- ⬜ RSM — no external reference checked yet

## Planned / wanted (from README, not started)
- ⬜ Statistics & ML: bootstrap CIs on fit parameters, outlier rejection on peak picks
- ⬜ Image Processing: 2D detector images (GIWAXS / RSM area images)
- ⬜ MATLAB Compiler standalone distribution (see Packaging above)

---
_Add new "wanted" items below as they come up so this stays the single source of truth._
