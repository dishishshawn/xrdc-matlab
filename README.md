# XRDC — X-ray Diffraction Analysis (Paik group)

MATLAB toolkit for X-ray diffraction analysis: reads Rigaku SmartLab and
PANalytical files, runs the common analyses (θ-2θ, rocking curve, φ scan, XRR,
reciprocal-space map), and produces publication-ready figures. A MATLAB port of
Dr. Tassilo Heeg's Delphi XRDC tool.

## Download & Install

### No MATLAB? (recommended for most users)

1. Go to the [latest release](https://github.com/dishishshawn/xrdc-matlab/releases/latest).
2. Download **`XRDCInstaller.exe`**.
3. Double-click it. On first run it downloads the matching **MATLAB Runtime**
   (a one-time, internet-required step — the Runtime is free, no MATLAB license
   needed).
4. Launch **XRDC** from the Start menu. The title bar shows the version, e.g.
   `XRDC Scan Analyzer v1.1.0`.

Windows only. The installer needs internet once to fetch the Runtime.

### Have MATLAB (R2022b or newer)?

1. Download the source (the **Source code (zip)** on the
   [latest release](https://github.com/dishishshawn/xrdc-matlab/releases/latest),
   or `git clone https://github.com/dishishshawn/xrdc-matlab.git`, or
   **Code → Download ZIP** for the latest `main`).
2. In MATLAB:
   ```matlab
   cd xrdc-matlab
   addpath(pwd)
   runtests      % optional — confirm green
   xrdcApp       % launch the GUI
   ```

Optional toolboxes (Signal Processing, Optimization, Curve Fitting, …) activate
better code paths automatically; every path has a pure-MATLAB fallback, so the
toolkit runs without them. See the toolbox table in
[`xrdc-matlab/README.md`](xrdc-matlab/README.md).

## Documentation

Full usage, scripts, the scan-struct schema, and the API are in
[`xrdc-matlab/README.md`](xrdc-matlab/README.md) and `xrdc-matlab/docs/`.

## Credits & citation

If you use this in a paper, cite Schwaigert et al. *J. Vac. Sci. Technol. A*
**41**, 022703 (2023) — the figure style is matched to that publication. The
underlying algorithms come from Dr. Tassilo Heeg's Delphi XRDC tool (FZJ/ISG1-IT).

Questions / bug reports: Shawn Agarwal (shawnagarwal0@gmail.com), Paik group.

## License

Academic use — see [`LICENSE`](LICENSE).
