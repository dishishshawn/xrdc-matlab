# xrdc-matlab — Scientific Assumptions for Deep Review

Each section below is a self-contained research prompt: paste into a literature-review agent or a domain expert without needing the repo. Goal is *rigorous external verification* before this toolkit is used to generate numbers for publication.

Tiers reflect the impact of being wrong:

- **Tier 1 — Numerical bias.** Results that go into a paper would be wrong.
- **Tier 2 — Convention drift.** Numbers may be self-consistent but mis-interpreted relative to other groups' published values.
- **Tier 3 — Code-level numerical correctness.** Small biases or edge-case bugs.

---

## TIER 1 — Could bias publication results

### 1.1 Kiessig fringe thickness ignores refraction (XRR)

**File:** `+xrdc/+lattice/thicknessFromFringes.m`

**Current implementation.** For N fringes at 2θ₁ < … < 2θ_N, the code computes both
```
t = (N − 1) · λ / [2 · (sin θ_N − sin θ_1)]            (closed form)
sin θ_i = m·i + b,  with t = λ / (2m)                  (linear fit)
```
i.e. it fits `sin θ` linearly in fringe index.

**The assumption.** This linear-in-`sin θ` formula is the **high-angle limit** of the Kiessig relation. It assumes refraction is negligible — that the X-ray refractive index inside the film equals unity. In real XRR, refraction makes the substrate-film system act as if the X-rays enter at a different angle than the goniometer reads. The correct relation is:
```
sin² θ_m − sin² θ_c = (m·λ / 2t)²
```
where θ_c is the film's critical angle for total external reflection. The correct fitting variable is `sin² θ` vs `m²`, not `sin θ` vs `m`.

**Research prompt.**

> The MATLAB function `thicknessFromFringes` fits XRR Kiessig fringes by linearising `sin θ_i` versus fringe index `i` and reporting thickness as `t = λ/(2·slope)`. For Cu Kα radiation on a typical perovskite oxide thin film (50–500 nm thick, fringes starting at 2θ ≈ 0.5°–2°), how large is the bias introduced by neglecting the critical-angle / refraction correction `sin² θ_m − sin² θ_c = (m·λ/2t)²` compared to the naive `sin θ_m ∝ m`? Provide:
> 1. The standard refraction-corrected Kiessig formula with citation (e.g. Tolan, *X-Ray Scattering from Soft-Matter Thin Films*; or Daillant & Gibaud).
> 2. A worked numerical example showing the percent error in reported thickness for a 100 nm PbTiO₃/SrTiO₃ film over the typical fringe-fit window (2θ ≈ 0.5°–3°).
> 3. Recommended replacement: the linear regression of `sin² θ_m` on `m²`, returning both `t` and `θ_c` simultaneously. Specify the fit weighting if fringes have non-uniform position uncertainty.
> 4. Whether published XRR papers (e.g. Schwaigert et al., *J. Vac. Sci. Technol. A* 41, 022703 (2023)) routinely apply the refraction correction, or whether the uncorrected formula is the lab convention.

---

### 1.2 Peak-fit uncertainties assume Gaussian, homoscedastic errors on count data

**File:** `+xrdc/+peaks/fitPeak.m` (lines 174–184)

**Current implementation.** `lsqcurvefit` minimises `Σ (y_i − f(x_i))²` (unweighted OLS). Parameter standard errors come from `σ² · (JᵀJ)⁻¹` with `σ² = RSS / (n − k)`.

**The assumption.** This formula is correct when residuals are **Gaussian and homoscedastic** (same variance everywhere). XRD counts are **Poisson** — variance equals the mean, so a 100,000-count peak top has σ ≈ 316 while a 100-count tail has σ ≈ 10. Unweighted OLS therefore:
- Over-weights peak tails relative to peak top
- Underestimates the SE on the centre and amplitude
- Produces R² values that aren't meaningful when intensities span orders of magnitude

The correct approach for Poisson counts is **weighted least squares** with `w_i = 1/y_i` (or `1/max(y_i, 1)` to handle zero bins), or equivalently maximum-likelihood Poisson fitting (Cash statistic).

**Research prompt.**

> A MATLAB XRD peak-fitter uses unweighted `lsqcurvefit` on raw count data, then computes parameter SE via the Jacobian as `σ² · (JᵀJ)⁻¹` with `σ² = RSS/(n−k)`. Counts are Poisson-distributed.
> 1. Quantify the bias on (a) the fitted peak centre 2θ, (b) the FWHM, (c) the amplitude, and their reported standard errors, for a Lorentzian on linear background where the peak amplitude is 10⁴ counts and the background is 10² counts.
> 2. Specify the correct weighted least-squares formulation: `w_i = 1/max(y_i, 1)` vs `w_i = 1/(y_i + 1)` vs the C-statistic / Cash statistic. Recommend one with citation (e.g., Bevington & Robinson; or Cash 1979 *ApJ* 228).
> 3. For the reported parameter SEs to be valid under WLS, what is the correct closed form? (Specifically: does `lsqcurvefit` with the `'Weights'` option propagate properly into Jacobian-based covariance, or does the implementation in MATLAB still need `σ² = 1` to interpret SEs correctly?)
> 4. Recommend whether existing Delphi XRDC fits in the literature would have to be re-run for parity, or whether the legacy unweighted approach is universally tolerated in XRD practice.

---

### 1.3 Hardcoded Cu Kα1 wavelength for all Rigaku files

**File:** `+xrdc/+io/readRigakuTxt.m:102`

**Current implementation.**
```matlab
scan.lambda = 1.5406;   % Cu Kα1, Rigaku SmartLab default (override if needed)
```
Every Rigaku `.txt` file gets λ = 1.5406 Å regardless of the optical configuration used.

**The assumption.** Rigaku ASCII exports do not contain the wavelength. The code assumes pure Cu Kα1, i.e. **Ge(220)×2 monochromator** on the SmartLab. But the same hardware operated in:
- **Bragg-Brentano** mode with no monochromator gives **Kα weighted average** ≈ 1.5418 Å (Kα1+Kα2, ratio 2:1)
- **HRXRD** with Ge(220)×2 gives essentially pure **Kα1** ≈ 1.5406 Å
- **HRXRD** with Ge(220)×4 also gives Kα1 but with better resolution (no wavelength effect)
- **PB (parallel-beam)** with parabolic mirror gives Kα average

Using 1.5406 Å for a Kα-mixed scan systematically shifts every reported d-spacing and lattice parameter.

**Research prompt.**

> Rigaku SmartLab ASCII (.txt) exports do not include the wavelength used. A MATLAB reader currently hardcodes λ = 1.5406 Å (Cu Kα1) for every file. Document:
> 1. The standard wavelengths in use today for Cu-source XRD: Kα1, Kα2, Kα weighted-mean (and the exact weighted ratio), and Kβ. Cite IUCr or NIST values to ≥ 5 significant figures.
> 2. Which Rigaku SmartLab optical configurations isolate pure Kα1 (Ge(220)×2 vs ×4 monochromator on the incident or diffracted side) versus configurations that pass the full Kα doublet. Confirm using Rigaku's own SmartLab manual or published method papers.
> 3. The propagated error on lattice parameter when using λ = 1.5406 Å where the data is actually Kα-weighted (≈ 1.5418 Å): for a SrTiO₃ (002) peak around 2θ ≈ 46.5°, what is the magnitude of the d-spacing offset and the implied lattice-parameter shift in Å?
> 4. Recommend whether a heuristic from filename or instrument metadata can disambiguate the configuration, or whether users should be required to supply `'Lambda'` explicitly when the configuration is not the lab default.

---

### 1.4 Pre-2019 physical constants in energy-to-wavelength conversion

**File:** `+xrdc/+lattice/energyToLambda.m` (lines 24–26)

**Current implementation.**
```matlab
h = 6.626068e-34;   % J·s    (pre-2019 CODATA, truncated)
c = 299792458;       % m/s   (exact)
e = 1.602e-19;       % C     (truncated to 4 sig figs)
```

**The assumption.** These are the values from the original Delphi source (matches `xrdc3.pas`). Since the 2019 SI redefinition:
- `h = 6.62607015 × 10⁻³⁴ J·s` (exact)
- `e = 1.602176634 × 10⁻¹⁹ C` (exact)
- `c = 299792458 m/s` (unchanged, exact)

This gives `hc/e = 12398.41984...` eV·Å (the code's docstring quotes `12398.4197`).

**Research prompt.**

> A MATLAB function uses the constants `h = 6.626068e-34 J·s`, `c = 299792458 m/s`, `e = 1.602e-19 C` to compute X-ray wavelength from photon energy. After the 2019 SI redefinition, the exact values are `h = 6.62607015e-34`, `e = 1.602176634e-19`. Compute:
> 1. The absolute and relative error in λ for a Cu Kα1 photon (~8047.8 eV) using the truncated vs SI-exact constants.
> 2. Whether this error propagates measurably into a typical XRD lattice-parameter measurement (uncertainty in `a` ≈ what fraction of a typical 0.001 Å measurement budget?).
> 3. Whether modern crystallography software (FullProf, GSAS-II, TOPAS) uses the SI-exact constants or carries the historical 12398.4 eV·Å value. Are there papers that document the conversion factor used in the calculation chain?
> 4. Recommend whether updating to SI-exact constants is worthwhile or whether the change is below measurement noise. Show the working.

---

### 1.5 Default Cu Kα1 line energy disagrees with NIST by 1.4 eV

**File:** `+xrdc/+data/xrayLines.json`

**Current implementation.**
```json
"Cu-Kalpha1-default": 8049.19,
"Cu-Kalpha1-NIST":    8047.8227,
```
The default value is selected by `"default": "Cu-Kalpha1-default"`, which differs from NIST by ~1.4 eV. The comment attributes this to "the historical value Dr. Heeg's program used."

**Research prompt.**

> A JSON of X-ray line energies includes both `Cu-Kalpha1-default = 8049.19 eV` and `Cu-Kalpha1-NIST = 8047.8227 eV` and selects the former as the default for all calculations.
> 1. What is the authoritative source for `8049.19 eV`? Is it an obsolete reference value (e.g. Bearden 1967), a Rigaku tube specification, or a misprint? Trace the provenance if possible.
> 2. What wavelength shift in Å does the 1.4 eV difference correspond to for Cu Kα1?
> 3. For a typical SrTiO₃ (002) peak at 2θ ≈ 46.5°, what is the implied d-spacing shift and the resulting lattice-parameter offset?
> 4. Which value is in use across the modern crystallography ecosystem (NIST X-ray Transition Energies Database; IUCr tables; FullProf default)?
> 5. Recommend whether the default should be switched to the NIST value, with a backwards-compatibility note about historical values used in this lab's prior publications.

---

### 1.6 RSM reciprocal-space convention: 2/λ versus 4π/λ

**File:** `+xrdc/+rsm/toReciprocalSpace.m` (line 107)

**Current implementation.**
```matlab
scale = 2 / lambda;
kPerp =  scale .* sin(θ) .* cos(ω − θ);
kPar  =  scale .* sin(θ) .* sin(ω − θ);
```

This uses **k = 2 sin θ / λ** (no factor of 2π) — the "Cullity convention" common in older XRD textbooks and in industry RSMs.

**The assumption.** Many modern publications (especially condensed-matter and thin-film growth papers) use **q = 4π sin θ / λ**, i.e. with 2π. The axes on a Schwaigert et al. JVST A 2023 RSM should be checked: are the k_par and k_perp axes in **Å⁻¹ (no 2π)** or **Å⁻¹ (with 2π)**?

**Research prompt.**

> A MATLAB reciprocal-space-map routine computes `k = (2/λ)·sin θ` for both in-plane and out-of-plane axes (no 2π factor — the Cullity convention).
> 1. Confirm the axis convention used in Schwaigert et al., *J. Vac. Sci. Technol. A* 41, 022703 (2023), Fig. 2(e). Is the colorbar axis `q_x, q_z (Å⁻¹)` with `q = 4π sin θ / λ`, or `k_x, k_z (Å⁻¹)` with `k = 2 sin θ / λ`? Verify against the paper's figure caption and methods section.
> 2. Tabulate which conventions are used by major thin-film MBE / XRD groups (e.g. Schlom group, Mannhart, Rijnders) and by RSM software packages (PANalytical Epitaxy/HRXRD, X'Pert Data Collector, Rigaku 3D Explore).
> 3. If the published target convention is `q = 4π sin θ / λ`, the current code's RSM peak positions are off by a factor of 2π — i.e. an SrTiO₃ (002) peak would land at the wrong place on the colorbar. Confirm this is a labeling issue (axis text wrong) rather than a numerical error (peak relative positions still correct), or vice versa.
> 4. Recommend whether to switch the default to `4π/λ` for consistency with modern papers, while keeping `2/λ` as an option for legacy work.

---

## TIER 2 — Convention drift / interpretation hazards

### 2.1 Nelson–Riley function mixes radians and degrees

**File:** `+xrdc/+lattice/nelsonRiley.m:50`

**Current implementation.**
```matlab
nrX = cos(θ).^2 ./ sin(θ) + cos(θ).^2 ./ thetaDeg;
```
where `θ` is in radians but `thetaDeg` is in degrees. The docstring calls this a "historical convention — do not 'fix' it."

**The assumption.** The classical Nelson–Riley function is `NR(θ) = ½·cos²θ·(1/sin θ + 1/θ)` with θ **in radians for both terms**. Mixing units leaves the **intercept** unchanged (since NR → 0 only when `cos θ → 0`, i.e. θ → 90°, regardless of which units θ is in for the dropping term) but **changes the slope**. The standard error of the intercept is therefore unaffected only if the units are consistent in the regression internals.

**Research prompt.**

> A Nelson–Riley extrapolation for lattice parameter uses the x-axis function `NR(θ) = cos²θ/sin θ + cos²θ/θ_deg` where the first `θ` is in radians but the second is in degrees. The intercept extrapolation point (θ = 90°) is independent of this unit mixing, but the slope and standard errors are not.
> 1. Verify the classical Nelson–Riley function: cite the original Nelson & Riley (1945) *Proc. Phys. Soc.* and a modern derivation (e.g. Cullity & Stock, *Elements of X-Ray Diffraction*, 3rd ed., Eq. 11-7). Confirm whether the standard form is `cos²θ·(1/sinθ + 1/θ)/2`, `cos²θ·(cos²θ/sin θ)` (the Taylor–Sinclair simplification), or the mixed form used here.
> 2. Quantify, for a typical lattice-parameter extrapolation with peaks spanning 2θ = 30°–120°, the difference in the *reported standard error* on `a_0` when using all-radians vs the current mixed-units form.
> 3. The code's intercept-SE formula uses the textbook OLS expression `Var(b_0) = σ²·Σxᵢ²/(n·Sxx)`, with `σ² = RSS/(n−2)`, deliberately departing from the Delphi original `xrdc3.pas:281`. Verify this is the correct simple-OLS formula and that no covariance correction is needed because the regression x-values are noise-free (only the y-values, the per-peak `a_i`, have uncertainty).
> 4. Recommend an unambiguous form: which form do modern crystallography textbooks (e.g. Pecharsky & Zavalij) use for the NR function today?

---

### 2.2 Pseudo-Voigt with shared FWHM is the Wertheim form, not TCH

**File:** `+xrdc/+peaks/fitPeak.m` (lines 305–317)

**Current implementation.**
```matlab
yL = lorentzian(x; x₀, fwhm)        % unit amplitude, same fwhm
yG = gaussian(x;  x₀, fwhm)         % unit amplitude, same fwhm
y = amp · (η·yL + (1 − η)·yG) + background
```

**The assumption.** This is the **Wertheim (1974)** form of the pseudo-Voigt, where both shapes share a single FWHM parameter. The fitted FWHM is a phenomenological width that doesn't separate the L and G contributions.

The modern standard for XRD is the **Thompson–Cox–Hastings (1987)** form, where the L and G components have *separate* FWHMs (Γ_L and Γ_G) that combine via a closed expression to give a **Voigt FWHM**, and η is a **derived** quantity:
```
Γ = (Γ_G⁵ + AΓ_G⁴Γ_L + BΓ_G³Γ_L² + CΓ_G²Γ_L³ + DΓ_GΓ_L⁴ + Γ_L⁵)^(1/5)
η = 1.36603(Γ_L/Γ) − 0.47719(Γ_L/Γ)² + 0.11116(Γ_L/Γ)³
```
This separates instrumental (Gaussian) and sample (Lorentzian) contributions cleanly.

**Research prompt.**

> A MATLAB XRD peak fitter implements pseudo-Voigt as `y = amp·(η·L + (1−η)·G)` with L and G sharing one FWHM parameter (Wertheim 1974).
> 1. Compare the Wertheim form to the Thompson–Cox–Hastings (1987) form, in which L and G have separate FWHMs combined into a Voigt FWHM. Cite both original papers.
> 2. For a thin-film diffraction peak where the instrumental broadening (Cu Kα1, Ge(220)×4 monochromator) is ~0.008° Gaussian and the sample broadening from crystallite size is ~0.05° Lorentzian, what difference does it make whether you report the fitted FWHM and η from Wertheim vs the TCH form?
> 3. Which form do FullProf, GSAS-II, and TOPAS use as the default profile in XRD Rietveld refinement, and which is dominant in single-peak fitting code?
> 4. Is the current Wertheim implementation safe for FWHM reporting in publications, or would the lab need to disclose this choice? Cite published examples either way.

---

### 2.3 adjustPeaks uses area bisector, not centroid

**File:** `+xrdc/+peaks/adjustPeaks.m` (lines 97–117)

**Current implementation.** The peak position is reset to the index that bisects the *area* of the half-height window — i.e. the cumulative-sum midpoint. The docstring explicitly notes this is "the Delphi convention, see ALGORITHM_SPEC §4.4" and "not the centroid, not the parabola vertex."

**The assumption.** For perfectly symmetric peaks, area-bisector and centroid coincide. For:
- **Asymmetric peaks** (e.g. unresolved Kα1/Kα2 doublet)
- **Peaks on a slope** (strong background gradient under the peak)
- **Peaks with shoulders or shoulders from adjacent reflections**

the area bisector and the *true* peak centre can differ by an appreciable fraction of the FWHM. The centroid (first moment) or a Gaussian/Lorentzian fit gives a position that is well-defined in terms of the underlying physical model.

**Research prompt.**

> A peak-position refinement routine replaces the initial peak 2θ estimate with the position where the cumulative area under the half-height window is split 50/50 (an area bisector, not a centroid).
> 1. Document the difference between the area bisector, the area centroid (first moment), the parabola vertex, and the peak maximum, for: (a) a pure Lorentzian, (b) a Cu Kα1+Kα2 unresolved doublet (separation Δλ/λ ≈ 0.0025, ratio 2:1), and (c) a Lorentzian on a linear background slope. Provide the relative shift in units of FWHM for each.
> 2. Which method is recommended for the most accurate lattice-parameter determination in the absence of full profile fitting? Cite Pecharsky & Zavalij or the IUCr Commission on Powder Diffraction newsletter.
> 3. Confirm whether the area-bisector method is used by any modern XRD software, or whether it is unique to the legacy Delphi XRDC. Is it documented in the published literature?
> 4. Recommend whether the centroid (`Σ x·y / Σ y` over the half-height window) is a strictly better default, and what backward-compatibility flag should be exposed to reproduce the Delphi behaviour.

---

### 2.4 simulatePattern does not apply structure factors

**File:** `+xrdc/+lattice/simulatePattern.m`

**Current implementation.** Enumerates `(h,k,l)` over the requested range, computes d-spacing, applies Bragg, and outputs every reflection that fits in the angular window. No structure-factor calculation; "forbidden" reflections (e.g. `(100)` for an FCC crystal, `(111)` for diamond) are *not* removed.

**The assumption.** For overlaying substrate ticks on a θ-2θ scan, the user is expected to know which (hkl) are allowed and to filter by index — e.g. requesting `[0 0 0 0 1 4]` to get only the (00ℓ) family of a tetragonal substrate. For a general crystal with a non-trivial basis, this would over-predict.

**Research prompt.**

> A MATLAB Bragg-position simulator computes peak 2θ for every (h, k, l) in a user-supplied range, ignoring structure-factor extinctions. It is used to overlay substrate ticks on θ-2θ scans of thin-film samples.
> 1. For the substrates currently in the lab's workflow — SrTiO₃ (Pm-3m, cubic), KTaO₃ (Pm-3m, cubic), LaAlO₃ (R-3c, rhombohedral pseudo-cubic), MgO (Fm-3m, cubic), Si (Fd-3m, diamond) — list the structure-factor-forbidden reflections that would currently be predicted as visible peaks if a user asked for the full (-3:3, -3:3, -3:3) range.
> 2. Quantify how often this matters in practice: for a SrTiO₃ (00ℓ) overlay (the common case), are there forbidden reflections within the 2θ range a researcher would typically display?
> 3. Recommend the minimal structure-factor logic that would catch the worst offenders without requiring a full structural model: e.g., systematic absences based on the lattice centering (P / F / I / C / R), not the basis. Cite the International Tables for Crystallography Vol. A, Ch. 2, for the absence rules.
> 4. For each substrate, list any (hkl) that would naively pass the simulator but be absent due to glide planes or screw axes in the space group, and the resulting visible (hkl)-set after filtering.

---

### 2.5 Background subtraction uses moving mean — peak-forest bias

**File:** `+xrdc/+signal/subtractBackground.m`

**Current implementation.** Default method is `'movmean'`. A moving average of *all* points — including the peaks themselves — is subtracted point-wise from the trace, and negative residuals are clipped to zero.

**The assumption.** This is fine when peaks are sparse and well-separated. For polycrystalline samples, near-edge regions, or scans with many overlapping reflections, the moving mean is **pulled up by the peaks** and the resulting baseline under-counts the true peak heights and depresses the area between peaks.

**Research prompt.**

> A MATLAB XRD background routine uses a moving-mean window of N points, subtracted pointwise from the counts, with negative results clipped to zero.
> 1. For a multi-peak XRD pattern (typical θ-2θ from a polycrystalline thin film: 5–10 peaks over 20°–80°), compare the recovered peak amplitudes after movmean-subtract vs the same trace processed with:
>    - Rolling minimum (e.g. `movmin`)
>    - Asymmetric Least Squares (AsLS — Eilers & Boelens 2005)
>    - Statistics-sensitive Non-linear Iterative Peak-clipping (SNIP — Morháč et al. *NIM A* 2000)
>    - The polynomial baseline used by FullProf and GSAS-II
>    Report bias on amplitude as a function of `windowSize / median(FWHM)` ratio.
> 2. Confirm that clipping negative residuals to zero (as the current code does) is standard practice or whether it biases reported integrated intensities, particularly when peaks are weak.
> 3. Recommend a default replacement that requires no extra user knob and works on both XRR (smooth long-range background) and θ-2θ (multi-peak) scans, or recommend separate defaults per scan type.

---

## TIER 3 — Code-level precision

### 3.1 Savitzky-Golay edge handling: docstring vs implementation mismatch

**File:** `+xrdc/+signal/derivatives.m` (lines 22–29, 66–67)

**Current implementation.**
```matlab
slope  = conv(counts, factorial(1) * flipud(g(:, 2)), 'same') / step;
slope2 = conv(counts, factorial(2) * flipud(g(:, 3)), 'same') / step^2;
```
But the docstring claims:
> the first and last (frameSize-1)/2 points of slope and slope2 use MATLAB's default sgolayfilt end behaviour, which is a lower-order polynomial fit to the available points.

`conv(..., 'same')` uses **zero-padding** at the boundaries, not a lower-order polynomial fit. `sgolayfilt` does a separate end-point polynomial fit. These produce different derivatives near the edges of the scan.

**Research prompt.**

> A MATLAB Savitzky-Golay derivative routine implements smoothing-derivatives via `conv(y, g, 'same')` but the docstring describes the end behaviour as if it used `sgolayfilt` (which does a lower-order polynomial extrapolation).
> 1. Quantify the difference at the first and last `(frameSize−1)/2` points between `conv` with zero-padding and `sgolayfilt` with end-point polynomial fits, on a typical XRD trace whose ends are not zero.
> 2. Recommend the appropriate edge handling for XRD-specific use cases: peak detection (where edge artifacts can spawn spurious peaks at the scan boundaries) vs. derivative-based background estimation.
> 3. Confirm whether the code should switch to `sgolayfilt(y, polyOrder, frameSize)` for the filtered signal, plus the matching derivative filters, to match the docstring claim. Specify the exact MATLAB call.

---

### 3.2 d-spacing formula verification: monoclinic, triclinic, rhombohedral

**File:** `+xrdc/+lattice/dSpacingFromHKL.m` (lines 85–108, 74–83)

**Current implementations.** All three formulas are derived; the rhombohedral and triclinic are reproduced inline above.

```
Monoclinic (β unique):
  1/d² = h²/(a·sinβ)² + k²/b² + l²/(c·sinβ)² − 2·h·l·cosβ/(a·c·sin²β)

Rhombohedral (a, α form):
  d² = a²·(1 − 3cos²α + 2cos³α) / [(h²+k²+l²)sin²α + 2(hk+kl+lh)(cos²α − cosα)]

Triclinic (a, b, c, α, β, γ):
  V²/(a²b²c²) = 1 − cos²α − cos²β − cos²γ + 2·cosα·cosβ·cosγ
  d² = V²/(a²b²c²) / [
        (h·sinα/a)² + (k·sinβ/b)² + (l·sinγ/c)²
      + 2·(h·k/(a·b))·(cosα·cosβ − cosγ)
      + 2·(k·l/(b·c))·(cosβ·cosγ − cosα)
      + 2·(h·l/(a·c))·(cosα·cosγ − cosβ) ]
```

**Research prompt.**

> Verify three d-spacing formulas in a MATLAB crystallography library against International Tables for Crystallography, Volume C, Table 1.2.1 (or equivalent authoritative source).
> 1. **Monoclinic, β-unique:** confirm the cross-term sign for `−2hl·cosβ/(a·c·sin²β)`. State the convention assumed for which angle is unique (β-unique is most common but some Russian/older European literature uses α-unique or γ-unique).
> 2. **Rhombohedral primitive setting (a, α):** verify the formula reproduced above. Confirm that the denominator group `2(hk + kl + lh)` has the correct order of indices and no sign change. Note: in some texts this is given as `2(hk + lh + kl)` with the same value, but if the formula uses `(cos²α − cosα)` instead of `(cosα − 1)cosα`, sign conventions can swap.
> 3. **Triclinic, general:** verify the cross-term signs `(cosα·cosβ − cosγ)`, `(cosβ·cosγ − cosα)`, `(cosα·cosγ − cosβ)`. These are commonly written with a + or − convention depending on whether the reciprocal-lattice metric tensor is computed via `g* = g⁻¹` or via direct cofactor expansion.
> 4. Cross-check by computing d for a known structure with all three formulas: e.g. quartz (trigonal-R or hexagonal setting), feldspar (triclinic, P-1, well-tabulated cell), and a monoclinic mineral such as gypsum (CaSO₄·2H₂O). For each test reflection, report the d computed by the formula and the published d. Flag any disagreement greater than 0.001 Å.

---

### 3.3 Auto prominence: log-domain criterion + unit-aware Poisson guard (replaced the 1.5% default, 2026-06-12)

**File:** `+xrdc/+peaks/findPeaks.m` (`autoDetect`, `estimateQuantum`)

**Current implementation.** When `MinProminence` is omitted, a peak must
(a) rise ≥ 0.3 decades above its neighbouring troughs on `log10(max(counts,1))`
(scale-invariant shape test) and (b) pass a Poisson significance test in *photon*
units: with `q` = the estimated count quantum, `N = counts/q`, `Nbg = bg/q`,
`bg = movmedian(counts, ~1°)`, require `z = (N − Nbg)/sqrt(N + Nbg + 1) ≥ 5`.
Explicit numeric `MinProminence` keeps the fixed linear threshold.

**The assumptions.**
- 0.3 decades (~2×) reflects how peaks are judged on the standard log-intensity
  plot. Calibrated on the Paik-group PTO/STO and PTO/LAO θ-2θ scans (the
  `validation/DATA_SWEEP.md` data dump) via `validation/runAutoProminence.m`.
- The guard is **unit-aware via the count quantum** `q` (`estimateQuantum`): raw
  photon-count data has `q = 1`; Rigaku counts-per-second exports show one photon
  as a small fixed intensity (~4.5–15 cps here) with a baseline dominated by exact
  zeros. When ≥ 5% of samples are exactly 0, `q` is the 1st percentile of the
  positive values (the single-photon floor); otherwise `q = 1`. This makes the 5σ
  test correct for both raw counts and cps — the earlier `5·sqrt(max(bg,1))` guard
  assumed raw counts and floored σ at 1 *count*, so on cps data (median background
  0) it became a flat 5-count bar that single-photon blips cleared, flooding real
  scans with hundreds of false peaks (e.g. 647 on S05).
- Using the peak's own `N` in the variance makes the test a difference-of-Poissons,
  slightly stricter than `sqrt(Nbg)` for strong peaks (which pass by orders of
  magnitude anyway) and decisive against single-photon blips (z ≈ 1).
- **Conservative by design:** the 5σ bar omits marginal ~3× background features.
  On weak-film scans (e.g. S05/S04/S10) this returns the substrate series only;
  genuinely strong films (S11, ~10× background) are detected. Users dig into weak
  films with an explicit `MinProminence`. This is a deliberate false-positive /
  false-negative trade chosen by the project owner.
- **Known fragility:** the "≥ 5% exact zeros" quantum trigger is a hard threshold.
  Data hovering near 5% zeros can flip `q` discontinuously (1 ↔ ~5) and shift which
  marginal peaks pass. Real Rigaku scans sit at 57–70% zeros — far from the cliff —
  so this is latent, not active; revisit (e.g. soft-blend or hysteresis) if a scan
  near the boundary ever misbehaves.

**Old default for reference (removed):** `max(1, 0.015·(max−min))` — missed
low-count film peaks next to 10⁶-count substrate peaks and over-segmented sharp
intense peaks (DATA_SWEEP findings, 2026).

---

### 3.4 XRR slab-model fitting (Parratt + Névot–Croce), 2026-06-13

**Files:** `+xrdc/+xrr/{opticalConstants,reflectivityModel,fitReflectivity}.m`,
`+xrdc/+data/atomicScattering.json`.

**Method.** Specular XRR is modelled by the Parratt recursion with Névot–Croce
interface roughness `r' = r·exp(−2 k_{z,j} k_{z,j+1} σ²)`. Optical constants come
from `δ,β = (rₑλ²/2π)·Σ Nⱼ f₁ⱼ,f₂ⱼ` using an embedded Henke/CXRO atomic-scattering
table **fixed at Cu Kα (8047.8 eV)** — other energies raise
`xrdc:xrr:unsupportedEnergy`. The fit minimises a **log-space** residual (XRR spans
~6 decades), seeded from `analyzeFringes` (thickness) and `findCriticalEdge`
(density), over film thickness, density, roughness, substrate roughness, scale,
background, and a footprint fill angle.

**Optimiser.** A multi-start bounded Nelder–Mead search over a fan of thickness
seeds (the fringe-spacing estimate plus a fixed spread), polished by `lsqnonlin`
(LM/trust-region); a Nelder–Mead-only path runs when the Optimization Toolbox is
absent. The multi-start is necessary, not decorative: single-seed LM gets trapped
in wrong fringe-count basins (recovering ~1 nm or ~55 nm instead of the true
thickness). Parameter SEs come from the Jacobian covariance via `pinv`, so the
identifiable parameters still get finite SEs when the scale/footprint columns are
near-degenerate.

**Assumptions / caveats.**
- Single film on a semi-infinite substrate; sharp slabs with Gaussian-blurred
  interfaces (Névot–Croce). No graded/SLD-profile layers, no multilayer (v1).
- Density is fit at the film's nominal composition (δ ∝ density); composition is
  not separable from density in XRR.
- The footprint correction `R·min(1, sinθ/sinθ_fill)` is a simple knife-edge
  illumination model; if the data is already footprint-corrected the fitted fill
  angle goes to ~0.
- **XRR is a phase problem: fits are non-unique.** The multi-start mitigates the
  thickness multimodality, but degenerate density/roughness/scale trade-offs
  remain. The fringe-spacing thickness (`analyzeFringes`) is an independent
  cross-check; large disagreement means distrust the fit.
- f₁/f₂ are tabulated only at Cu Kα; anomalous-dispersion fine structure near
  absorption edges is not modelled beyond those fixed values.

---

## Summary

| Tier | Item | Lines |
|------|-----|--------|
| 1 | Kiessig refraction correction | thicknessFromFringes.m |
| 1 | Poisson errors / unweighted SE | fitPeak.m |
| 1 | Hardcoded Rigaku λ = 1.5406 | readRigakuTxt.m:102 |
| 1 | Pre-2019 constants | energyToLambda.m:24-26 |
| 1 | Cu Kα1 default 8049.19 eV | xrayLines.json |
| 1 | RSM 2/λ vs 4π/λ | toReciprocalSpace.m:107 |
| 2 | Nelson–Riley mixed units | nelsonRiley.m:50 |
| 2 | Pseudo-Voigt Wertheim vs TCH | fitPeak.m:305-317 |
| 2 | adjustPeaks area-bisector | adjustPeaks.m:97-117 |
| 2 | No structure-factor extinctions | simulatePattern.m |
| 2 | movmean background bias | subtractBackground.m |
| 3 | Savitzky-Golay edge handling | derivatives.m:66-67 |
| 3 | d-spacing formula verification | dSpacingFromHKL.m:85-108 |
| 3 | Log-domain auto prominence + unit-aware Poisson guard | findPeaks.m (autoDetect) |
| 3 | XRR slab-model fitting (Parratt + Névot–Croce) | fitReflectivity.m, reflectivityModel.m, opticalConstants.m |

Before peer distribution, **Tier 1.3 (Rigaku wavelength)** and **Tier 1.5 (Cu Kα1 default 8049.19 vs NIST 8047.8)** are the highest priority because they affect every Rigaku file users will throw at the tool with no obvious diagnostic.
