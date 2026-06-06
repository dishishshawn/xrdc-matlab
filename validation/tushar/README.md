# Validation set — Tushar

Cross-check the xrdc-matlab app against a known-good result from Tushar: run his raw
scan through the app and confirm the output matches the graph he produced.

## Folders
| Folder        | What goes here                                                          |
| ------------- | ----------------------------------------------------------------------- |
| `input/`      | The raw scan file(s) Tushar sent (`.txt` Rigaku, `.xrdml`, `.x00`, …).   |
| `reference/`  | The graph/figure he sent that the data should reproduce (PNG/PDF/etc.).  |
| `output/`     | What you export from the app — compared against `reference/`.            |

Contents of these three folders are git-ignored (they're his data, not a deliverable);
the folder structure, this README, and `RESULTS.md` are tracked.

## How to run the comparison

1. Drop Tushar's scan file into `input/` and his graph into `reference/`.
2. In MATLAB, from the inner project folder (`xrdc-matlab/`):
   ```matlab
   addpath(pwd)
   xrdcApp
   ```
3. Click **Load Scan...**, pick the file from `input/`. The app auto-detects the scan
   type (θ-2θ / rocking curve / φ / XRR / RSM) and runs the matching analysis.
4. Match the on-screen settings to whatever Tushar used (axis scaling, fit shape,
   background, peak thresholds) so it's an apples-to-apples comparison.
5. Click **Export 600 dpi...** and save into `output/` with a name that mirrors the
   input (e.g. `S04_2theta_omega_xrdc.png`).
6. Open `output/` and `reference/` side by side and compare — see the checklist below.

Prefer scripting? Point a demo at the file instead of the GUI:
```matlab
fname = '../../validation/tushar/input/<his-file>';   % path is relative to xrdc-matlab/
run('examples/demoThetaTwoTheta.m')                   % or the demo matching the scan type
```

## What to compare
- [ ] Peak **positions** (2θ / ω) line up within instrument resolution.
- [ ] Peak **intensities / relative heights** match.
- [ ] Fit-derived numbers agree: FWHM (rocking curve), film thickness in nm (XRR
      Kiessig fringes), d-spacing / lattice parameter (θ-2θ).
- [ ] Axis ranges, log vs linear Y, and any substrate overlay match his figure.
- [ ] Note any deliberate divergences (see `xrdc-matlab/docs/USER_GUIDE.md` §7) before
      calling a difference a bug.

Record each run in `RESULTS.md` so the comparison is reproducible.
