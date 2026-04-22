---
title: XRDC Algorithm Specification
subtitle: Reference for the MATLAB port
source: Delphi/Pascal source in ./XRD Converter Source (T. Heeg, FZJ/ISG1-IT, 2003+)
status: draft
---

# XRDC Algorithm Specification

This document is a working reference for porting XRDC (X-Ray Diffraction analysis, originally written in Delphi by Tassilo Heeg) to MATLAB. It extracts the algorithms, equations, thresholds, and edge cases from the original source so the MATLAB implementation can reproduce the original results where that is desirable, and deliberately diverge where modern MATLAB tooling is clearly better.

Every equation below was taken directly from the Pascal source; line references are given in the form `filename:line` to make verification easy.

Scope of the port (Dr. Paik's goals, for context):

1. Direct read of Rigaku `.raw` binary files (not supported by the original).
2. Cleaner, more organised data handling in MATLAB.
3. Consistent, publication-ready figure formatting across datasets.

A pragmatic rule of thumb for every algorithm in this document: replicate the original's numerical result where it is already correct, modernise the implementation (e.g. closed-form OLS, `lsqcurvefit`, vectorised operations) where the original was constrained by Pascal's lack of numerical libraries.

---

## 1. Data model

The original program revolves around three types (`xrdc1.pas:16-66`):

`TDataPoint` — one diffraction point: `twotheta`, `counts`, `slope`, `slope2` (first and second derivative in 2θ).

`TDataSet` — one scan. Holds:

- `OrigDataSet` — untouched copy of the raw data, for reset.
- `DataSet` — current data, possibly after smoothing / background subtraction.
- `DataSize` — number of points.
- `Identifier` — short label for the legend.
- `SecondAxis` / `SecondAxisName` — the constant axis value for θ-2θ and similar scans (e.g. the Omega value at which a θ-2θ was taken).
- `PeakSet` — list of `TPeak` records.

`TScanType` — `stTwoTheta_Omega`, `stPhi`, `stOmega`, `stPsi`, `stTheta` (constant-ω θ-2θ). The scan type drives file parsing, axis labelling, and RSM eligibility.

MATLAB mapping. A single `struct array` `scans(i)` with fields `{twoTheta, counts, slope, slope2, identifier, secondAxis, secondAxisName, scanType, peaks, original}` is sufficient. Keep an immutable `.original` so that noise suppression and background subtraction can be undone non-destructively — the original Delphi program relies on this behaviour.

---

## 2. File-format parsers

`ReadFile` (`xrdc1.pas:438`) is a dispatcher on the first non-empty line of the file:

- Starts with `XML VERSION` or `<XRDMEASUREMENTS` → PANalytical XML (`.xrdml`).
- Starts with `HR-XRDSCAN` → Philips legacy text (`.x00`, High-Resolution format).
- Starts with `SUPERDATA` → Picker (very old; prompts user to select a dataset from a list).
- Anything else → two-column plain text.

Important: the loader always treats a comma and a period as interchangeable decimal separators via `DotCommaToDecimalSep` / `DecimalSepToDot` (`xrdctool.pas:13-28`) because the original was developed under German Windows locale. **The MATLAB port should continue to accept both.**

### 2.1 PANalytical `.xrdml` (XML)

`ReadFileTypeXML` (`xrdc1.pas:494`) uses `TXMLParser` (from `LibXmlParser.pas`) with XPath-like traversal. The parser:

1. Scans for `scanAxis` attribute on `xrdMeasurements/xrdMeasurement/scan` and maps it to `TScanType` by string match: `"2Theta-Omega"`, `"Omega"`, `"Phi"`, `"Psi"`, `"Reciprocal Space"` (area scan with 2Theta-Omega slices).
2. Reads `startPosition` and `endPosition` for the primary axis; computes `stepwidth = (end - start) / (N - 1)` where N is the number of intensity values (derived from whitespace-separated counts).
3. Reads `SecondAxis` (constant axis) from the scan's auxiliary axis element.
4. For each data point `i`, `twoTheta(i) = start + i * stepwidth` and `counts(i) = parseInt(token_i)`.

Edge case: reciprocal-space area maps are ingested as multiple θ-2θ slices with different `SecondAxis` (ω-offset) values. These are the scans that feed the RSM viewer (§7).

MATLAB port note: Use `readstruct` with `"FileType","xml"` (R2020b+) or the older `xmlread` DOM API. The XML has a few vendor-specific namespaces; strip the namespace prefixes or use local-name matching.

### 2.2 Philips `.x00` (text)

`ReadFileTypePhilips` (`xrdc1.pas:991`) — line-based key/value header followed by a SCANDATA block:

- Keys of interest: `FIRSTANGLE`, `STEPWIDTH`, `NROFDATA`.
- Block marker: `SCANDATA` — subsequent lines are integer counts, one per line.
- Same `start + i*stepwidth` reconstruction of 2θ as above.

### 2.3 Picker `.1035`, `.596` (text, historical)

`ReadFileTypePicker` (`xrdc1.pas:1065`). Begins with `SUPERDATA`, contains multiple concatenated scans. Columns (right-to-left) are `H K L D TIM CTS TTH CHI PHI OM` — only `CTS` (counts) and `TTH` (2θ) are kept. The program prompts the user to select a scan by index. Rarely used; the MATLAB port can implement this last (or skip entirely) unless a specific dataset requires it.

### 2.4 Plain text (two-column)

`ReadFileTypeText` (`xrdc1.pas:1235`). Two columns: 2θ (float) and counts (float). Rules:

- If any count exceeds 2×10⁹, divide the entire column by 10 — a heuristic to undo fixed-point × 10 encoding sometimes seen in vendor exports.
- Compute `stepwidth = 2θ[1] - 2θ[0]`, then assert it is uniform within 0.001°. On failure, warn the user but continue (the original issues a dialog; the MATLAB port should log a warning).
- Counts are rounded to integers (the original storage type).

### 2.5 Rigaku `.raw` (NEW — not in original)

Not implemented in the Delphi source. The Rigaku `.raw` binary has two common variants: the legacy "RAW1.01" and the newer "RAS" ASCII. Suggested ingestion path for the MATLAB port:

- Detect variant from the magic: `"RAW1.01"` or `"RAW1.02"` as ASCII at offset 0 → binary; `"*RAS_DATA_START"` → ASCII.
- For the ASCII RAS variant, parse the `*MEAS_COND_*` key/value header (especially `MEAS_COND_AXIS_NAME`, `MEAS_SCAN_START`, `MEAS_SCAN_STOP`, `MEAS_SCAN_STEP`, `MEAS_SCAN_SPEED`, `MEAS_SCAN_AXIS_X`, `MEAS_SCAN_UNIT_X`, `MEAS_SCAN_UNIT_Y`) and the `*RAS_INT_START` / `*RAS_INT_END` data block.
- Map the scan axis string to `TScanType` in the same way as XRDML (§2.1).
- The binary RAW variant is layout-documented in Rigaku SmartLab application notes; a standalone reader function should live in `+xrdc/+io/readRigakuRaw.m`.

Target contract for every parser in the MATLAB port:

```
function scan = readScan(path)
% returns struct with fields:
%   twoTheta       (Nx1 double, degrees)
%   counts         (Nx1 double, raw counts — no rounding)
%   scanType       (string: "twoThetaOmega" | "omega" | "phi" | "psi" | "theta" | "area")
%   secondAxis     (scalar double, value of constant axis; NaN if not applicable)
%   secondAxisName (string)
%   identifier     (string, from filename)
%   metadata       (struct with full header for traceability)
```

---

## 3. Signal processing

### 3.1 Slopes (first and second derivative)

`CalcSlopes` (`xrdc1.pas:1387`). Sliding-window linear fit over ±`SpinEditSlopeCnt` points (default 5–10); the fitted slope is stored as `DataSet[i].slope`, and a second OLS pass over the slopes yields `slope2`. The helper is `Geradenanpassung` (`xrdc1.pas:1353`), a closed-form OLS:

```
Sx  = Σ xᵢ,  Sy  = Σ yᵢ
Sxx = Σ xᵢ²,  Sxy = Σ xᵢ yᵢ
m   = (n·Sxy − Sx·Sy) / (n·Sxx − Sx²)
b   = (Sxx·Sy − Sx·Sxy) / (n·Sxx − Sx²)
```

Port note: the Delphi code returns `(x0, y0)` where `x0` is the slope `m` and `y0` is the intercept `b`. Do not let the variable names fool you when cross-checking.

MATLAB: `movmean` is not appropriate here because it is a mean, not a slope. Use `movslope` (signal processing toolbox) or `conv` with a Savitzky–Golay derivative kernel. A single call to `sgolayfilt(x, order, frame)` with `deriv=1` replaces `CalcSlopes`; the same with `deriv=2` replaces the second pass and is numerically better behaved because it fits the derivative analytically.

### 3.2 Noise suppression (moving average)

`NoiseSuppression` (`xrdc1.pas:1420`). Symmetric moving-average smoothing of `counts` with window `m` (user-controlled). Endpoints use the available points only (narrow window near the edges). Not applied to `slope` / `slope2`; those are recomputed afterwards.

MATLAB: `smoothdata(counts, 'movmean', m)` or `movmean(counts, m, 'Endpoints','shrink')`.

### 3.3 Background subtraction

`SubstractBackground` (`xrdc1.pas:1457`). A moving-average baseline with a larger window (`SpinEdit5`, typically 50–200 points); the baseline is subtracted pointwise. Negative values are clipped to zero.

MATLAB: same idea. A documented alternative is `msbackadj` (Bioinformatics Toolbox) or a rolling-minimum / asymmetric-least-squares method; the simplest port keeps the original moving-average behaviour for reproducibility and exposes a second, optional "modern" method behind a `method` name–value argument.

---

## 4. Peak detection

Two detectors run on the smoothed-and-background-subtracted trace. Both are driven by the pre-computed `slope` / `slope2`.

### 4.1 ScanPeak1 — slope zero-crossing (`xrdc1.pas:1482`)

Walk left-to-right through the trace:

1. Find a local maximum in `slope` exceeding a threshold `T` (`SpinEdit6`, in counts / degree).
2. Find the zero-crossing `slope → 0`.
3. Confirm with a subsequent local minimum in `slope` with value `< -maxSlope / 2`.
4. Peak position is the average of the max-slope index and the min-slope index (not the zero-crossing index itself — this reduces bias from asymmetric peaks).

### 4.2 ScanPeak2 — second-derivative (`xrdc1.pas:1542`)

Peak candidate when `slope` changes sign **and** `slope2 < -SpinEdit6`. Equivalent to finding inflection-bounded concave-down regions.

### 4.3 FindPeaks (`xrdc1.pas:1566`)

Driver that runs either ScanPeak1 or ScanPeak2 (user choice, `Button14` calls with algorithm 2). Post-processing:

- Reject any candidate with `counts < MinPeak` (UI-controlled absolute threshold).
- Sort by 2θ ascending.
- Merge peaks closer than `FloatEdit13` (default 0.05°) — keep the higher-intensity one.

### 4.4 AdjustPeaks — FWHM-based refinement (`xrdc1.pas:1715`)

For every accepted peak, the position is refined using 100× interpolation of the raw data around the peak:

1. Interpolate `counts(2θ)` to a grid 100× finer than the native step.
2. Determine local background from the edges of the peak window.
3. Compute FWHM with a user-controlled threshold `SpinEdit3` (percent of height above background; default 50%).
4. Peak position is the midpoint of the FWHM-interval — this is the area bisector, not the centroid.
5. After refinement, collapse duplicates that are within 0.03° of each other (post-refinement pairs sometimes collide).

MATLAB port guidance: `interp1(..., 'pchip')` for step 1; `findpeaks` with `'MinPeakProminence'` and `'WidthReference','halfheight'` does most of 3–4 natively. The MATLAB port can keep the ScanPeak1/2/AdjustPeaks pipeline for faithful reproduction, but the default user-facing path should be `findpeaks` with sensible prominence defaults because it is more robust and exposes uncertainty directly.

### 4.5 Substrate-corrected peak adjustment (Nelson–Riley feedback loop)

`AdjustPeaksBySubstrateNelsonRiley` (`xrdc1.pas:2909`) re-calibrates peaks against a known substrate after a Nelson–Riley fit (§6.2):

```
a_corr = y0 + x0·(cos²θ / sinθ  +  cos²θ / θ_deg)
nλ     = 2·a_corr·sinθ
θ_new  = arcsin(nλ / (2·y0))
```

where `(x0, y0)` are the slope and intercept of the Nelson–Riley line. The mixed `cos²θ/sinθ + cos²θ/θ_deg` term is Nelson–Riley's empirical extrapolation (note: the `θ_deg` in the denominator is *in degrees*, not radians — historical convention, see `xrdc3.pas`; the MATLAB port should comment this loudly to avoid silent unit errors).

---

## 5. Peak fitting (Lorentz / Gauss)

`TForm8.DoFit` in `xrdc9.pas` implements a brute-force 3-parameter grid search for a single peak:

**Model.** Background is linear `bg(x) = bga·x + bgb`, fitted from averages of the trace in the edge regions `bgposleft` and `bgposright` (positions to the left/right of the peak, user-set):

```
y1 = mean(counts in window around (bgposleft + 2θmin)/2)
y2 = mean(counts in window around (bgposright + 2θmax)/2)
bga = (y2 - y1) / (bgposright - bgposleft)
bgb = y1 - bga · bgposleft
```

**Peak shape (user choice).**

- Lorentz:  `f(x) = yscale / (1 + (xscale·(x−x0))²)`  →  `FWHM = 2 / xscale`.
- Gauss:    `f(x) = yscale · exp(−(xscale·(x−x0))²)`  →  `FWHM = 2·√(−ln 0.5) / xscale`.

**Search.** Over `x0 ∈ [2θmin, 2θmax]`, `xscale ∈ [1/(2θmax−2θmin), 10/step]`, `yscale ∈ [(peak−bg)/3, 2·peak]`:

- 20 grid points per parameter on each pass (20³ = 8000 evaluations).
- Fit metric: sum of squared residuals (implicit).
- Refine the region around the best triple and repeat.
- Stop when `|devmin − devbest| < devbest / 1e10` **and** the `x0` step is below 0.001°.

MATLAB port guidance: replace the entire search with `lsqcurvefit` (Optimization Toolbox) or `fit` with `'lorentzian1'` / `'gaussian1'` (Curve Fitting Toolbox). Seed initial values from the brute-force grid's first pass or from `findpeaks` output. Preserve the legacy brute-force path behind `method = "bruteforce"` only if bit-for-bit reproduction of old analyses is needed.

Return values from the fit (for consistency with the original UI):

- `x0` (peak 2θ, degrees)
- `FWHM` (degrees)
- `yscale` (peak height above background)
- `bga`, `bgb` (background slope/intercept)
- Residual standard deviation

---

## 6. Lattice-parameter calculation

All of Section 6 is in `xrdc3.pas`.

### 6.1 Physical constants and Bragg's law

```
h  = 6.626068e-34 J·s
c  = 299792458 m/s
e  = 1.602e-19 C
λ[Å] = 1e10 · (c·h) / (E[eV] · e)       (wavelength from energy)
d    = λ / (2 · sin θ)                  (Bragg's law, n=1 implied)
```

`CalculateHKL` computes the multiplier from an expected lattice parameter `a_expected` the user typed:

```
m_i = round(a_expected / d_i)
a_i = d_i · m_i
```

### 6.2 Simple arithmetic mean

`CalculateArithLattice`: `a_mean = mean(a_i)`; standard error via `CalculateStandardDeviation`:

```
σ = √( Σ(aᵢ − a_mean)² / ((n−1)·n) )
```

(This is the standard error of the mean, *not* the sample std-dev — worth preserving the formula but renaming the reporting in the UI to `SE(a)`.)

### 6.3 Nelson–Riley extrapolation

`CalculateNelsonRiley`. For each peak, compute the NR x-axis value:

```
NR(θ) = cos²θ / sinθ  +  cos²θ / θ_deg         (θ_deg is θ in degrees)
```

Linear fit `a_i  =  y0  +  x0 · NR(θ_i)` using the same closed-form OLS as §3.1. The intercept `y0` is the extrapolated lattice parameter at θ = 90° (where NR → 0), which is the best estimate free of most systematic errors (specimen displacement, absorption, zero offset).

Regression standard errors (textbook OLS):

```
SE(y0) = √( Sxx · RSS / ((n−2)·(n·Sxx − Sx²)) )
SE(x0) = √(  n  · RSS / ((n−2)·(n·Sxx − Sx²)) )
```

where `RSS` is the residual sum of squares, `Sx = Σ NR(θᵢ)`, `Sxx = Σ NR(θᵢ)²`, `x0` is the slope and `y0` the intercept.

**Delphi discrepancy:** `xrdc3.pas:281` writes the intercept SE with an extra `n` in the denominator — `SE_delphi(y0) = √(Sxx·RSS / (n·(n−2)·(n·Sxx − Sx²)))` — i.e. it underestimates the SE by a factor of √n. This is a latent bug in the Delphi source; the MATLAB port uses the textbook formula and does **not** reproduce the Delphi scaling (per CLAUDE.md "match the algorithm, not the implementation"). The slope SE has no Delphi counterpart (only `deltay0` is computed there).

MATLAB implementation: `polyfit` + `polyval` gives the fit; `fitlm` gives the uncertainty directly. Either works; `fitlm` is preferred for the output table.

### 6.4 Energy / line selection

`XRAY.def` enumerates the X-ray lines:

```
Cu Kα1 default   8049.19  eV       (the standard for perovskite work)
Cu Kα1 NIST      8047.8227
Cu Kα2 NIST      8027.8416
Cu Kα           8038.5           (weighted α1/α2 average)
Cu Kβ1 NIST      8905.413
W  Lα1 NIST      8398.242        (present for tungsten-anode tubes)
...
```

The wavelength computed from these energies is the `λ` threaded through Bragg's law and the RSM transform.

---

## 7. Reciprocal-space maps

### 7.1 Transform (`Prepare2dPlot`, `xrdc1.pas:3223`)

For each point in every θ-2θ slice with its own `SecondAxis` (ω). The transform is slightly asymmetric in how it applies the 2θ offset:

```
θ_raw    = 2θ_point / 2                              (uncorrected, for building ω)
ω        = SecondAxis − (2θ_center / 2) + θ_raw + ΔΩ_RSM        (in radians)
θ        = (2θ_point + ΔΘ_RSM) / 2                   (corrected θ, for Bragg)
k_perp   = (2/λ) · sin(θ) · cos(ω − θ)
k_par    = ±(2/λ) · sin(θ) · sin(ω − θ)              (sign flips with `RSMPlotFlip`)
```

Units: both `k` components in Å⁻¹. `2θ_center` is the mean of the first and last 2θ of the slice, used to normalise the ω baseline so different slices share a consistent k-space origin.

Subtlety preserved from the Delphi source (`xrdc1.pas:3293–3296`): `θ` is used in two different senses. The `θ_raw` that builds `ω` is the *uncorrected* half-angle, but the `θ` in `sin(θ)·cos(ω−θ)` is the 2θ-corrected one. Reproduce this asymmetry — undoing it silently will shift RSM peaks off-target when users have non-zero `ΔΘ_RSM`.

For non-θ-2θ scans (e.g. ω-scans at fixed 2θ), the simpler branch is used:

```
ω = (SecondAxis + ΔΩ_RSM) · π/180
θ = (2θ_point + ΔΘ_RSM) · π/360
```

### 7.2 Offsets (`SetRSMOffsets`, `xrdc1.pas:4020`)

The user clicks on a known peak (typically the substrate's symmetric reflection) in the 2D viewer. The program:

1. Finds the slice whose `SecondAxis` is closest to the click's y-coordinate.
2. Computes `domega = SecondAxis − (2θ_center / 2)` for that slice.
3. Opens a dialog (Form12) pre-filled with `BaseTwoTheta = x_click`, `BaseOmega = x_click/2 + domega`, `BaseOmegaOffset = domega`.
4. The user enters the known theoretical 2θ/ω for that peak; the dialog computes and returns `ΔΘ_RSM`, `ΔΩ_RSM`, which the transform then applies to every point.

Port note: this offset mechanism is essential for publication-quality RSMs. Without it, the absolute positions are off by the goniometer's zero error. Keep the interactive "click-to-offset" UX in the MATLAB port — it is by far the most user-friendly way to do this.

### 7.3 2D rendering

The original has two 2D renderers:

- `TwoDPlot.pas` — hand-rolled raster to a `TCanvas`, with three colour modes (`0` colour, `1` B&W, `2` inverted B&W). Used on-screen.
- Gnuplot export via `xrdc_gnuplot_template.plt` — PNG at 2560×1920, contour-on-surface, custom palette `rcol(x), gcol(x), bcol(x)` (a piecewise linear RGB colour map defined in the template).

The template has `#PNGONLY#` / `#SCREENONLY#` / `#COLORONLY#` / `#BWONLY#` / `#NOCONTOURONLY#` tokens that are string-substituted in `xrdc11.pas` before the file is passed to `wgnuplot.exe`.

MATLAB port: collapse both paths into a single `imagesc` / `contourf` call with a named colour map. For publication output:

- Consistent figure size, font, colour map across all scans (this is goal #3 of the project).
- The original's colour palette is essentially a blue → cyan → green → yellow → red ramp; MATLAB's `turbo` is a close analogue, or define a matching palette by implementing `rcol/gcol/bcol` from lines 17–19 of the gnuplot template and verifying visually on a reference RSM.

---

## 8. Structure simulation (Bragg peak ticks)

`CalculateStructureLines` (`xrdc1.pas:4341`). Given lattice parameters `a,b,c,α,β,γ` and an hkl range, generate predicted 2θ positions for the selected crystal system:

- **Cubic** (`ComboBox3=0`):  `d = a / √(h² + k² + l²)`
- **Tetragonal** (`1`):       `1/d² = (h² + k²)/a² + l²/c²`
- **Orthorhombic** (`2`):     `1/d² = h²/a² + k²/b² + l²/c²`
- **Hexagonal** (`3`):        `1/d² = (4/3)·(h²+k²+hk)/a² + l²/c²`    (reported as 4-index `h (-(h+k)) k l`)
- **Rhombohedral** (`4`):     `d² = a²·(1 − 3cos²α + 2cos³α) / [(h²+k²+l²)·sin²α + 2(hk+kl+lh)(cos²α − cosα)]`
- **Monoclinic** (`5`):       `1/d² = h²/(a sin β)² + k²/b² + l²/(c sin β)² − 2hl·cos β / (a·c·sin²β)`
- **Triclinic** (`6`):        the full general form (see source for numerator/denominator with `cos α cos β cos γ`).

Then Bragg: `2θ = (360/π) · arcsin(λ/(2d))` (degrees).

Post-processing:

1. Bubble-sort by `d` descending (equivalent: by `2θ` ascending).
2. Merge duplicates within `|Δd| < 1e-6 Å`. Tie-break by first `min(h+k+l)`, then `min(l·10000 + k·100 + h)`.
3. Write results to `StrucSimPeaks[]` with `(h, k, l, d, 2θ, description)`.

MATLAB port: vectorise over the hkl grid (no triple loop). Crystal-system dispatch stays a `switch` on a system-enum. The "merge duplicates" rule is there because symmetry-equivalent reflections (e.g. {100}/{010}/{001} in cubic) share a `d`; you may prefer to leave them in and rely on multiplicity instead, but matching the original output means emitting only one.

Substrate-reference patterns (`Substrates.def`, German-locale-encoded) are pre-computed 2θ values for a set of common substrates (SrTiO₃, LaAlO₃, MgO, DyScO₃, NdGaO₃, Si, CeO₂, GaN, BaTiO₃, Au, GdScO₃, SrO) and include an `ALLOWED=` line that selects which orders to display (e.g. `1,2,3,4` → only 100, 200, 300, 400 for the `(x00)` families). Port these as a shipped resource file (`+xrdc/+data/substrates.json` or similar).

---

## 9. Film thickness from Laue / Kiessig fringes

`EstimateThickness` (`xrdc15.pas`) uses the selected fringe peaks and fits `d` from the spacing:

```
d[nm] = (n − 1) · λ[nm] / (2 · (sin θ_n − sin θ_1))
```

where `θ_1 … θ_n` are the selected fringe peaks in ascending 2θ, `n` is the number of peaks, and `λ` is in nm (the code divides the Å wavelength by 10). This is the Laue equation rearranged for period-from-N-fringes.

A parallel `CalculateHKL` in the thickness form supports half-integer fringe indices (checkbox 2), useful for some fringe-numbering conventions where the zeroth fringe is between the first two observed peaks:

```
m_i = baseindex + i            (integer mode)
m_i = baseindex + i + 0.5      (half-integer mode)
a_i = d_peak,i · m_i            (d from Bragg for each fringe)
```

MATLAB port: direct translation. Use `polyfit` of `sin(θ_i)` against `i` to get thickness from the slope — equivalent, and gives a proper uncertainty estimate for free.

---

## 10. Plot formatting (goal #3)

The original's on-screen plot uses the custom `graphplot.pas` GDI renderer; PNG export uses gnuplot with the template in §7.3. For the MATLAB port, a single `+xrdc/+plot/plotScan.m` function should produce publication-ready figures with these invariants:

- **Fonts**: Arial / Helvetica, sizes 18 (ticks) / 20 (axis labels) / 22 (title) — matching the gnuplot template's screen sizes.
- **Line width**: 1.5 for main trace, 1.0 for additional scans, 1.0 dashed for structure ticks.
- **Background**: white; grid off by default; minor ticks on.
- **Axis labels**: `2\theta [\circ]` for X; `Counts` (log scale by default for XRD) for Y.
- **Log-Y default**: the original clamps `counts ≤ 0` to `1` before display (`PaintGraph`, `xrdc1.pas:1901-1904`). Keep this behaviour to avoid `log(0)` in the MATLAB port.
- **Multi-scan stacking**: the original multiplies each additional scan by `3^(j+1)` (or, for >50 scans, by `10^(20·(j+1)/N)`) to vertically separate them on a log plot (`xrdc1.pas:1939-1953`). This is a purely cosmetic "waterfall" offset — preserve the behaviour but expose the factor as a parameter.
- **Colour**: deterministic palette for the first N scans (not random). The original actually generates random RGB per additional scan (`Round(Random($FFFFFF))`) which is not publication-friendly; replace with a deterministic palette (e.g. `lines` or a colour-blind-safe sequence like `turbo(N)` / `parula(N)` sampled at even intervals).
- **Export**: `exportgraphics(fig, 'out.png', 'Resolution', 600)` for raster; `exportgraphics(fig, 'out.pdf', 'ContentType','vector')` for vector.

---

## 11. Things to drop from the original

- **Auto-updater** (`Check4SoftwareUpdate`, HTTP GET of `currentversion.xml`). Not relevant to MATLAB.
- **Gnuplot dependency** (`xrdc_gnuplot_template.plt`, `xrdc.INI → [GNUPLOT]`). MATLAB's native graphics replaces it.
- **Random per-scan colours**. Replaced with deterministic palette above.
- **German-locale decimal-separator rewriting** (`DotCommaToDecimalSep`) — still needed on read, but internal storage is always `double` with `.` as the separator. Only the writer paths (if any CSV export) need to worry.
- **Forms for modal UI** (Delphi-specific). MATLAB's port can expose the same workflows through `uifigure`/App Designer or a plain scripting API. For Dr. Paik's goals, a scripting API with a thin App Designer wrapper is likely better than reproducing the 17-form dialog tree.

---

## 12. MATLAB package layout (recommendation)

```
+xrdc/
  +io/
    readScan.m           % dispatcher
    readXrdml.m
    readPhilipsX00.m
    readRigakuRaw.m      % NEW
    readRigakuRas.m      % NEW (ASCII variant)
    readTextScan.m
  +signal/
    smooth.m             % movmean wrapper
    subtractBackground.m
    slopes.m             % Savitzky-Golay derivatives
  +peaks/
    findPeaks.m          % modern default, uses built-in findpeaks
    findPeaksLegacy.m    % ScanPeak1/2 for reproduction
    adjustPeaks.m        % FWHM-based refinement
    fitPeak.m            % lsqcurvefit with lorentzian/gaussian choice
  +lattice/
    bragg.m              % d <-> 2theta <-> lambda
    nelsonRiley.m
    dSpacingFromHKL.m    % all 7 crystal systems
  +rsm/
    toReciprocalSpace.m  % the Prepare2dPlot transform
    setOffsetsInteractive.m
  +plot/
    plotScan.m
    plotStack.m          % waterfall
    plotRsm.m
    publicationStyle.m   % applies the invariants in §10
  +data/
    substrates.json      % from Substrates.def
    xrayLines.json       % from XRAY.def
```

---

## 13. Verification plan

Before declaring the MATLAB port correct, run each of the following against the Delphi binary's output on the same input:

1. Parse one PANalytical `.xrdml` file; compare `twoTheta` and `counts` vectors element-for-element.
2. Smooth with the same window; compare to 1e-12.
3. Run peak detection with identical thresholds; expect the same peak list (FWHM refinement may differ at the 5th decimal because of interpolation method).
4. Nelson–Riley fit on a 5-peak Si calibration set; compare intercept and slope to 1e-6.
5. Peak-fit a test Lorentzian with known FWHM; compare fitted FWHM to 1e-3 (both programs should nail a perfect Lorentzian).
6. Generate a structure-simulation pattern for cubic SrTiO₃ with λ=1.5406 Å, hkl ∈ [−4,4]; compare 2θ list.
7. RSM transform a small test area scan; compare `(k_par, k_perp)` for every point to 1e-6 Å⁻¹.

Then, for Rigaku-specific validation (no Delphi equivalent):

8. Ingest three `.raw` files the lab has on hand; visually check that 2θ range, step size, and total counts match what Rigaku's own software reports.

---

*End of spec. This document will be revised as the MATLAB port progresses; see commit history for changes.*
