# xrdc-matlab

MATLAB toolkit for X-ray diffraction analysis in the Paik group. Reads Rigaku SmartLab and PANalytical files, runs the common analyses (θ-2θ, rocking curve, φ scan, XRR, reciprocal-space map), and produces publication-ready figures matching Schwaigert et al. *J. Vac. Sci. Technol. A* **41**, 022703 (2023), Fig 2.

Two ways to use it: a **GUI** for point-and-click work, or **scripts** when you want full control.

---

## 1. Install

You need **MATLAB R2022b or newer**.

**Recommended toolboxes.** Every code path has a pure-MATLAB fallback, so the toolkit still runs without any of these. Installing them activates the better path automatically.

| Toolbox                              | Activates                                                                                       |
| ------------------------------------ | ----------------------------------------------------------------------------------------------- |
| Signal Processing Toolbox            | `findpeaks` (XRD/Kiessig peak detection), `sgolay`/`sgolayfilt` derivatives, XRR log-envelope   |
| Optimization Toolbox                 | `lsqcurvefit` peak fits with Jacobian-based parameter SEs                                       |
| Curve Fitting Toolbox                | `csaps` smoothing-spline option in `smoothCounts` and `subtractBackground`; `fit` + `confint` for Kiessig thickness CIs |
| Global Optimization Toolbox          | `Method="multistart"` in `fitPeak` — robust against poor initial guesses                        |
| Parallel Computing Toolbox           | `parfor` over files in `loadAreaScan` (RSM I/O speed-up)                                        |
| Statistics and Machine Learning      | Bootstrap CIs on fit parameters, outlier rejection on peak picks (planned)                      |
| Image Processing Toolbox             | Only needed for 2D detector images (GIWAXS/RSM area scans) — not yet used                       |
| MATLAB Compiler                      | Only needed if distributing the GUI as a standalone executable                                  |

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

`xrdc.io.readScan` auto-detects from the first 512 bytes (or the extension for binary formats):

- **Rigaku SmartLab** `.txt` (both headered and headerless variants)
- **Rigaku SmartLab Studio II** `.hgx` (GlobalFit project — an HDF5 container; reads the measured curve, plus the simulated fit into `metadata`)
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

**License:** academic use — see [`LICENSE`](../LICENSE).

---

## 7. Standalone executable (run without a MATLAB license)

The GUI can be compiled into a double-clickable **`XRDC.exe`** so collaborators who don't have a MATLAB license can run it. This needs **MATLAB Compiler** on the *build* machine only — end users install just the free MATLAB Runtime.

**Build it** (from the repo root, on a machine with MATLAB Compiler + the optional toolboxes installed):

```matlab
addpath build
buildStandalone                 % → build/standalone/XRDC.exe
buildStandalone(Embed=true)     % also bundles the Runtime installer (larger)
```

**Run it** (end-user machine, no MATLAB license required):

1. Install the **MATLAB Runtime** — the version must match the MATLAB release used to build (`mcrversion` on the build box). Free download: <https://www.mathworks.com/products/compiler/matlab-runtime.html>. This is a multi-GB, one-time install.
2. Double-click `XRDC.exe`. The GUI opens exactly as `xrdcApp` does in MATLAB.

**Notes & caveats:**

- The Runtime is **license-free but not small** — going standalone removes the *license* requirement, not the disk footprint.
- The compiler bakes in whichever code path the build machine resolves, so **build on a fully-licensed box** to get the Signal Processing / Optimization / Curve Fitting paths rather than the pure-MATLAB fallbacks.
- The exe is **Windows-only**. For macOS/Linux, build on that platform (see the note in `build/buildStandalone.m`).
- Each build is tied to one Runtime version; ship the matching Runtime (or use `Embed=true`).
- The MATLAB *project* remains the source of truth — the exe is a release artifact, not a replacement.
