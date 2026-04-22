# xrdc-matlab

MATLAB port of XRDC — the X-ray diffraction analysis tool originally written in Delphi by Dr. Tassilo Heeg (FZJ/ISG1-IT, 2003+). Reimplemented for the Paik group to handle Rigaku SmartLab exports, organise data more cleanly, and produce consistent publication-quality figures matching Schwaigert et al. *J. Vac. Sci. Technol. A* 41, 022703 (2023).

See `../ALGORITHM_SPEC.md` for algorithm derivations, `../PROJECT_PLAN.md` for the phased roadmap, and `docs/USER_GUIDE.md` for usage.

## Status

**Phases 1–6 complete. Handoff-ready pending MATLAB run of the full test suite on the lab's machine.**

Implemented:

- `+xrdc/+io/` — PANalytical XRDML (XML), Philips `.x00`, plain text, **Rigaku SmartLab `.txt`** (both headered and headerless variants). Format-sniffing dispatcher (`readScan`) with scanType inference from filename for Rigaku files.
- `+xrdc/+signal/` — moving-average smoothing, background subtraction, Savitzky-Golay 1st/2nd derivatives.
- `+xrdc/+peaks/` — prominence-based detector (`findPeaks`, wraps MATLAB `findpeaks`), legacy slope / slope² detectors (`findPeaksLegacy`, direct port of Delphi `ScanPeak1`/`ScanPeak2`), FWHM area-bisector refinement (`adjustPeaks`), Lorentz/Gauss/pseudo-Voigt fitting (`fitPeak`, `lsqcurvefit` by default with Jacobian SEs; Delphi 20³ grid search under `'Method','bruteforce'`).
- `+xrdc/+lattice/` — Bragg's law, energy↔wavelength, d-spacing for all seven crystal systems, Nelson-Riley extrapolation (textbook OLS SEs — deliberately diverges from `xrdc3.pas:281`, see ALGORITHM_SPEC §6.3), Kiessig-fringe film thickness, vectorised structure simulation with symmetry-equivalent merging.
- `+xrdc/+rsm/` — reciprocal-space coordinate transform (θ-asymmetry preserved per ALGORITHM_SPEC §7.1), folder/file-list area-scan loader, click-to-align goniometer-offset workflow.
- `+xrdc/+plot/` — `publicationStyle` (central style knob for the JVST-A-2023 invariants), `plotScan`, `plotStack` (waterfall with the Delphi 3^(j) / 10^(20·j/N) multiplier and a deterministic palette), `plotRsm` (filled contour with log decade colorbar — Fig 2(e) style).
- `+xrdc/+data/` — `xrayLines.json`, `substrates.json` (ported from `XRAY.def` / `Substrates.def`).
- Unit tests for every package under `tests/`. Tests requiring real lab data auto-skip via `assumeTrue(isfile(...))` when the data folder is absent.
- `examples/` — six demo scripts covering θ-2θ, XRR, rocking curve, phi scan, RSM, and structure-simulation workflows.
- `docs/USER_GUIDE.md` + `docs/RIGAKU_NOTES.md`.

Stubbed (not in the Paik lab's current workflow — kept behind clear `xrdc:io:notImplemented` error messages in case historical archive files surface):

- `xrdc.io.readRigakuRas` — ASCII RAS (`*MEAS_COND_*` header + `*RAS_INT_START` block)
- `xrdc.io.readRigakuRaw` — binary RAW1.01/RAW1.02

GUI: `xrdcApp.m` — self-contained `uifigure` app for lab members who prefer point-and-click. Auto-detects scan type from the file, runs the matching analysis, live plot preview, parameter tweaks re-run on the fly, one-click 600 dpi export. Invoke with `xrdcApp` after adding the repo to the path. See *Quick start* below.

## Requirements

- MATLAB R2022b or newer (uses `arguments` blocks and `smoothdata`).
- **Signal Processing Toolbox** — `sgolay`, `sgolayfilt`, `findpeaks`.
- **Optimization Toolbox** — `lsqcurvefit` in `xrdc.peaks.fitPeak`; fallback `'Method','bruteforce'` works without it.
- **Statistics and Machine Learning Toolbox** — recommended (used by `fitlm` in one Nelson-Riley cross-check test).

## Quick start

### GUI (for lab members)

```matlab
cd xrdc-matlab
addpath(pwd)
xrdcApp                                           % opens the interactive app
```

Click *Load Scan...*, pick any `.txt` or `.xrdml`. The app detects the scan type, runs the appropriate analysis, shows a live preview, and lets you tweak parameters on the fly. *Export 600 dpi...* writes a publication-ready PNG.

### Scripts (for power users)

```matlab
cd xrdc-matlab
addpath(pwd)
runtests                                          % full test suite
run('examples/demoThetaTwoTheta.m')               % θ-2θ demo
run('examples/demoRsmKTaO3.m')                    % RSM demo — JVST A 2023 Fig 2(e) style
```

All demos export PNG at 600 dpi into the current directory. To run a demo on your own file, set `fname` in the workspace first:

```matlab
fname = 'my_scan.txt';
demoRockingCurve                                  % uses your file
```

## Paper parity

Target style: Schwaigert et al. *J. Vac. Sci. Technol. A* 41, 022703 (2023), Fig 2.

- Fig 2(a) θ-2θ log-Y with substrate ticks — `demoThetaTwoTheta.m`
- Fig 2(b) zoom + Laue fringes — combine with `xlim` on the output of `plotScan`
- Fig 2(c)/(d) rocking curves + FWHM — `demoRockingCurve.m`
- Fig 2(e) RSM contourf + log decade colorbar — `demoRsmKTaO3.m`
- Fig S1 XRR + Kiessig thickness — `demoXRR.m`

## Package layout

```
xrdc-matlab/
├── +xrdc/
│   ├── +io/        file-format parsers + dispatcher (XRDML, Philips .x00, Rigaku .txt, plain text)
│   ├── +signal/    smoothing, background, Savitzky-Golay derivatives
│   ├── +peaks/     prominence + legacy detectors, FWHM refine, peak fit
│   ├── +lattice/   Bragg, d-spacing, Nelson-Riley, structure simulation, film thickness
│   ├── +rsm/       reciprocal-space transform, area-scan loader, click-to-align
│   ├── +plot/      publicationStyle, plotScan, plotStack, plotRsm
│   └── +data/      substrates.json, xrayLines.json
├── examples/       6 runnable demo scripts
├── tests/          matlab.unittest cases (per-package)
├── docs/           USER_GUIDE.md, RIGAKU_NOTES.md
├── xrdcApp.m       GUI (uifigure) — load → auto-analyze → export
├── runtests.m      test-suite entry point
├── README.md       this file
└── CLAUDE.md       conventions for future AI contributions
```

## Dropping features from the original

The following Delphi-era features are intentionally **not** ported — see `../PROJECT_PLAN.md` §2:

- Auto-updater (HTTP version check).
- Gnuplot PNG export (replaced by MATLAB's native `exportgraphics`).
- Obsolete Picker `.596` / `.1035` format (`readScan` raises `xrdc:io:notSupported`).
- German-locale-specific `.def` / `.INI` file writing (we read them, but store as JSON).
- Random per-scan colours (replaced with deterministic palette).

## References

- Original source: `../XRD Converter Source/` (Delphi).
- Algorithm spec: `../ALGORITHM_SPEC.md`.
- Project plan: `../PROJECT_PLAN.md`.
- Paper reference: Schwaigert et al. *J. Vac. Sci. Technol. A* 41, 022703 (2023).
