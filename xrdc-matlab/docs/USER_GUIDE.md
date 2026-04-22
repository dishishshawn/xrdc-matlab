# xrdc-matlab User Guide

MATLAB port of Dr. Tassilo Heeg's Delphi XRDC, built for the Paik lab. Reads Rigaku SmartLab ASCII (`.txt`), PANalytical XML (`.xrdml`), Philips `.x00`, and plain text; produces publication-ready figures matching the Paik group's paper style (Schwaigert et al. *J. Vac. Sci. Technol. A* 41, 022703 (2023)).

This guide covers installation, the data model, and the five core workflows:

1. Load a scan
2. Plot a single θ-2θ / XRR trace
3. Rocking curve → FWHM in arcsec
4. φ scan → substrate symmetry check
5. Reciprocal-space map (RSM) from multiple slices

For algorithm derivations see `../../ALGORITHM_SPEC.md`. For contribution conventions see `../CLAUDE.md`.

---

## 1. Installation

### Requirements

- **MATLAB R2022b or newer** (uses `arguments` blocks and `smoothdata`).
- **Signal Processing Toolbox** — required (`sgolay`, `sgolayfilt`, `findpeaks`).
- **Optimization Toolbox** — required for `xrdc.peaks.fitPeak` (`lsqcurvefit`). Legacy brute-force grid is available as a fallback via `'Method','bruteforce'`.
- **Statistics & Machine Learning Toolbox** — recommended (used by `fitlm` in the Nelson–Riley cross-check test).
- **Curve Fitting Toolbox** — optional; not currently required.

### Path setup

From MATLAB, inside the `xrdc-matlab/` directory:

```matlab
cd xrdc-matlab
addpath(pwd)             % makes the +xrdc package visible
runtests                  % full test suite — expect green
```

No other third-party dependencies. Nothing is compiled. The package is pure MATLAB code under `+xrdc/`.

### Repository layout

```
xrdc-matlab/
  +xrdc/
    +io/         — format parsers and scan-struct helpers
    +signal/     — smoothing, background subtraction, derivatives
    +peaks/      — peak detection and fitting (modern + legacy paths)
    +lattice/    — Bragg, d-spacing, Nelson-Riley, simulatePattern, thicknessFromFringes
    +rsm/        — reciprocal-space coordinate transform, area-scan loader
    +plot/       — plotScan, plotStack, plotRsm, publicationStyle
    +data/       — xrayLines.json, substrates.json
  examples/      — one script per workflow (run from repo root)
  tests/         — matlab.unittest suites (one per package)
  docs/          — this guide + Rigaku format notes
```

---

## 2. The scan struct

Every reader returns the same struct, defined by `xrdc.io.emptyScan()`:

| Field            | Type           | Meaning                                                |
| ---------------- | -------------- | ------------------------------------------------------ |
| `twoTheta`       | Nx1 double     | Primary axis (degrees). For ω-scans / RC this is ω.    |
| `counts`         | Nx1 double     | Intensity, raw counts (not rounded)                    |
| `scanType`       | string         | `"twoThetaOmega"`, `"omega"`, `"phi"`, `"psi"`, `"theta"`, `"area"`, `"unknown"` |
| `secondAxis`     | scalar double  | Constant axis value (e.g. ω for a θ-2θ slice). `NaN` if N/A |
| `secondAxisName` | string         | Name of the constant axis                              |
| `identifier`     | string         | Short label (filename stem by default)                 |
| `lambda`         | scalar double  | Wavelength in Å, or `NaN` if not stored                |
| `metadata`       | struct         | Full header for traceability                           |
| `sourcePath`     | string         | Original file path                                     |
| `sourceFormat`   | string         | `"xrdml" | "philipsX00" | "text" | "rigakuTxt"`        |

All public APIs take this struct. Angles in degrees; wavelengths and lattice parameters in Å; thicknesses in nm.

---

## 3. Loading scans

```matlab
scan = xrdc.io.readScan('path/to/file');
```

`readScan` sniffs the format from the first 512 bytes:

| Signature                             | Format            |
| ------------------------------------- | ----------------- |
| `<?xml` or `<xrdMeasurement`          | PANalytical XRDML |
| `HR-XRDSCAN`                          | Philips `.x00`    |
| `*RAS_DATA_START` or `*MEAS_COND`     | Rigaku RAS (stubbed — not in lab workflow) |
| `RAW1.0…`                             | Rigaku binary (stubbed) |
| `INTENSITY, CPS` anywhere in header   | Rigaku SmartLab `.txt` |
| Filename starts with `TR_`            | Rigaku SmartLab `.txt` (fallback cue) |
| anything else                         | plain-text two-column |

The dispatcher picks the right sub-parser; you call `readScan` for everything. Individual parsers (`readXrdml`, `readRigakuTxt`, etc.) can also be invoked directly.

---

## 4. Workflows

All five demo scripts live in `examples/`. Run from the repo root:

```matlab
run('examples/demoThetaTwoTheta.m')
```

### 4.1 θ-2θ scan with substrate overlay

`examples/demoThetaTwoTheta.m` — loads a PbTiO₃/SrTiO₃ θ-2θ, detects peaks, overlays simulated SrTiO₃ 00L positions.

Core calls:
```matlab
scan = xrdc.io.readScan(path);
pk   = xrdc.peaks.findPeaks(scan, 'MinProminence', 100);
T    = xrdc.lattice.simulatePattern( ...
          struct('system','cubic','a',3.905), [0 0 0 0 1 4], scan.lambda, ...
          'TwoThetaRange', [scan.twoTheta(1), scan.twoTheta(end)]);
h    = xrdc.plot.plotScan(scan, 'Title', "PbTiO_3/SrTiO_3 — θ-2θ");
```

### 4.2 XRR with Kiessig-fringe thickness

`examples/demoXRR.m` — specular reflectivity, fringe detection, thickness from periodicity.

```matlab
pk    = xrdc.peaks.findPeaks(subScan, 'MinProminence', 0.02*max(sub.counts));
thick = xrdc.lattice.thicknessFromFringes([pk.twoTheta], scan.lambda);
fprintf('d = %.1f ± %.1f nm\n', thick.thicknessFitNm, thick.thicknessFitSeNm);
```

`thicknessFromFringes` fits `sin(θ_i)` vs fringe index and returns both the N-fringe closed form and the fit thickness with uncertainty.

### 4.3 Rocking curve — FWHM in arcsec

`examples/demoRockingCurve.m` — Lorentzian fit on a ±0.5° window; converts FWHM from degrees to arcsec.

```matlab
pk     = xrdc.peaks.findPeaks(scan, 'MinProminence', 0.1*max(scan.counts));
fit    = xrdc.peaks.fitPeak(scan, [pk(1).twoTheta - 0.5, pk(1).twoTheta + 0.5], ...
            'Shape', 'lorentz');
arcsec = fit.fwhm * 3600;
```

Note: Rigaku exports label the RC x-axis `"2θ"` even though the data are ω. `readRigakuTxt` sets `scanType = "omega"` from the filename (`*RC*`, `*rocking*`); demo relabels the plot axis accordingly.

### 4.4 φ scan — 4-fold symmetry check

`examples/demoPhiScan.m` — linear-Y plot (log-Y hides off-peak noise), detects peaks spaced ≥30°, reports spacings.

### 4.5 Reciprocal-space map

`examples/demoRsmKTaO3.m` — reproduces the Fig 2(e) style of the KTaO₃ MBE paper: filled contour on log intensity, decade-tick colorbar (1, 10, 10², …, 10⁵).

```matlab
scans = xrdc.rsm.loadAreaScan('path/to/folder', 'Pattern', '*.xrdml');
h     = xrdc.plot.plotRsm(scans, ...
          'Mode',       "contourf", ...
          'Imin',       1, ...
          'Imax',       1e5, ...
          'Colormap',   "turbo", ...
          'ExportPath', 'out.png');
```

The transform in `xrdc.rsm.toReciprocalSpace` preserves the θ-asymmetry documented in ALGORITHM_SPEC §7.1 (`θ_raw` builds ω; corrected `θ` feeds the k formulas). Do not override this — it matters when applying goniometer zero-offset corrections.

For interactive goniometer-offset alignment, use `xrdc.rsm.setOffsetsInteractive` — click the known substrate peak in the RSM figure, enter the theoretical 2θ/ω, and the function returns `(ΔΘ, ΔΩ)` for use with `plotRsm`.

---

## 5. Plot customisation

All plotting functions funnel through `xrdc.plot.publicationStyle`. To change the look globally (font family, tick/label sizes, grid defaults), edit that function; callers don't need to change.

Per-plot name/value pairs on `plotScan` / `plotStack` / `plotRsm`:

- `FontName` — default `"Arial"`
- `TickFontSize` / `LabelFontSize` / `TitleFontSize` — default 18 / 20 / 22
- `LogY` — default `true` (clamps counts ≤ 0 to 1, mirroring `PaintGraph` in `xrdc1.pas:1901`)
- `LineWidth` — default 1.5 for the main trace
- `ExportPath` — pass a filename to auto-save via `exportgraphics` at 600 dpi

Paper-style defaults are locked to match Schwaigert et al. JVST A 41 (2023), Fig 2.

---

## 6. Testing

Full suite:

```matlab
results = runtests('tests');
assertSuccess(results)
```

Individual suite:

```matlab
runtests('tests/testRsm.m')
```

Tests that depend on real Paik-lab data files gate on file existence via `assumeTrue(isfile(...))` — they skip (not fail) when the data folder is absent. These are listed under each file's header.

Verification plan (ALGORITHM_SPEC §13): V1–V7 can run from tests/; V8 (Rigaku 2θ range / step / counts match Rigaku's own software) is verified by the `testReadRigakuTxtRealFiles` integration test.

---

## 7. Known divergences from Delphi XRDC

These are deliberate — documented in ALGORITHM_SPEC and enforced in tests:

| Topic                          | Delphi XRDC                        | xrdc-matlab                                 |
| ------------------------------ | ---------------------------------- | ------------------------------------------- |
| Nelson–Riley SE (intercept)    | `xrdc3.pas:281` formula (off by 1/√n) | Textbook OLS (see ALGORITHM_SPEC §6.3)  |
| Peak fitting default           | 20³ brute-force grid               | `lsqcurvefit` with Jacobian SEs. `'Method','bruteforce'` keeps the legacy path. |
| `plotStack` colours            | Random RGB per trace               | Deterministic palette                       |
| Auto-updater / gnuplot export  | Present                            | Dropped                                     |
| Picker `.596` / `.1035` files  | Supported                          | Dropped (obsolete)                          |
| German-locale decimal on write | Yes                                | Never write comma decimals                  |

---

## 8. Getting help

- Read the function docstring first (`help xrdc.peaks.fitPeak`).
- For algorithm questions, check the referenced line number in the original Delphi source.
- For format questions see `docs/RIGAKU_NOTES.md` (Rigaku) and ALGORITHM_SPEC §2 (all formats).
- Ownership: Shawn Agarwal (shawnagarwal0@gmail.com), Paik group stakeholder.
