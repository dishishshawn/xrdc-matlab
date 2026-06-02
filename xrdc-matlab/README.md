# xrdc-matlab

MATLAB toolkit for X-ray diffraction analysis in the Paik group. Reads Rigaku SmartLab and PANalytical files, runs the common analyses (θ-2θ, rocking curve, φ scan, XRR, reciprocal-space map), and produces publication-ready figures matching Schwaigert et al. *J. Vac. Sci. Technol. A* **41**, 022703 (2023), Fig 2.

Two ways to use it: a **GUI** for point-and-click work, or **scripts** when you want full control.

---

## 1. Install

You need **MATLAB R2022b or newer**.

**Recommended toolboxes** (code has fallbacks when missing, but performance and fit quality improve with them):

| Toolbox                              | Used for                                                              |
| ------------------------------------ | --------------------------------------------------------------------- |
| Curve Fitting Toolbox                | Peak fits, Parratt/fringe fits, smoothing splines for log-detrend     |
| Optimization Toolbox                 | Constrained `lsqnonlin`/`fmincon` for layer parameters                |
| Signal Processing Toolbox            | `findpeaks`, filtering, FFT-based fringe period extraction            |
| Global Optimization Toolbox          | `MultiStart`/`ga`/`particleswarm` for multilayer XRR fits             |
| Parallel Computing Toolbox           | `parfor` batch refinement, GPU Parratt recursion                      |
| Statistics and Machine Learning      | Bootstrap CIs on fit params, outlier rejection on peak picks          |
| Image Processing Toolbox             | Only if handling 2D detector images (GIWAXS/RSM area scans)           |
| MATLAB Compiler                      | Only if distributing the GUI as a standalone executable               |

Most university TAH licenses include all of these — check `ver` in MATLAB before requesting additions.

Clone or copy the repo, then from MATLAB:

```matlab
cd xrdc-matlab
addpath(pwd)
runtests                                          % optional — confirm green
```

That's it. Nothing is compiled; everything lives under `+xrdc/`.

---

## 2. Use the GUI (recommended if you don't write MATLAB)

```matlab
xrdcApp
```

1. Click **Load Scan...** and pick any `.txt` (Rigaku SmartLab) or `.xrdml` (PANalytical) file.
2. The app detects the scan type (θ-2θ, rocking curve, φ scan, XRR) and runs the matching analysis automatically.
3. Adjust parameters on the left panel — the plot re-runs live.
4. Click **Export 600 dpi...** to save a publication-ready PNG.

That covers most day-to-day work. For RSM and anything custom, drop down to the scripts.

---

## 3. Use the scripts

Each common workflow has a demo in `examples/`. Run a demo as-is to see the expected output, then point it at your own file by setting `fname` in the workspace first:

```matlab
fname = 'my_scan.txt';
run('examples/demoRockingCurve.m')
```

| Workflow                       | Demo script                  | What it does                                                  |
| ------------------------------ | ---------------------------- | ------------------------------------------------------------- |
| θ-2θ scan + substrate overlay  | `demoThetaTwoTheta.m`        | Loads scan, finds peaks, overlays simulated substrate 00L     |
| Rocking curve → FWHM           | `demoRockingCurve.m`         | Lorentzian fit, FWHM in arcsec                                |
| φ scan symmetry check          | `demoPhiScan.m`              | Linear-Y plot, detects 4-fold peaks, reports spacings         |
| XRR + Kiessig thickness        | `demoXRR.m`                  | Specular fit, fringe detection, film thickness in nm          |
| Reciprocal-space map           | `demoRsmKTaO3.m`             | Filled contour, log decade colorbar (Fig 2(e) style)          |
| Structure simulation           | `demoStructureSim.m`         | Predict 2θ for a given crystal system + (hkl) list            |
| Rigaku I/O sanity              | `demoRigakuWorkflow.m`       | End-to-end read → smooth → peaks → plot on a Rigaku file      |

All demos export 600 dpi PNGs to the current directory.

If you'd rather call the functions directly:

```matlab
scan = xrdc.io.readScan('path/to/file.txt');     % auto-detects format
pk   = xrdc.peaks.findPeaks(scan, 'MinProminence', 100);
fit  = xrdc.peaks.fitPeak(scan, [22.5 23.5], 'Shape', 'lorentz');
xrdc.plot.plotScan(scan, 'Title', "PbTiO_3/SrTiO_3");
```

**See `docs/USER_GUIDE.md` for the full API, the scan-struct schema, and per-workflow detail.**

---

## 4. Supported file formats

`xrdc.io.readScan` auto-detects from the first 512 bytes:

- **Rigaku SmartLab** `.txt` (both headered and headerless variants)
- **PANalytical** `.xrdml` (single scan and multi-scan area files)
- **Philips** `.x00`
- Plain two-column text

Rigaku binary `.raw` and ASCII `.ras` are not currently in the lab workflow and are stubbed — let me know if you have files in those formats and we'll wire them in.

---

## 5. Layout

```
xrdc-matlab/
  +xrdc/          core package — io, signal, peaks, lattice, rsm, plot, data
  examples/       7 demo scripts (run from repo root)
  tests/          matlab.unittest suites
  docs/           USER_GUIDE.md, RIGAKU_NOTES.md
  xrdcApp.m       GUI entry point
  runtests.m      test-suite entry point
```

---

## 6. Help and citation

- **Function docs:** `help xrdc.peaks.fitPeak` (or any other function) from the MATLAB prompt.
- **Workflow detail:** `docs/USER_GUIDE.md`.
- **Rigaku format notes:** `docs/RIGAKU_NOTES.md`.
- **Bug reports / questions:** Shawn Agarwal (shawnagarwal0@gmail.com), Paik group.

If you use this in a paper, please cite Schwaigert et al. *J. Vac. Sci. Technol. A* **41**, 022703 (2023) — the figure style is matched to that publication, and the original Delphi XRDC tool by Dr. Tassilo Heeg (FZJ/ISG1-IT) is the source of the underlying algorithms.
