# Design: XRR slab-model fitting (Parratt + Névot–Croce)

**Date:** 2026-06-13
**Status:** Approved (Shawn, via design review in session)
**Feature:** A full reflectivity forward model + fitter in `+xrdc/+xrr` that fits a
single film on a substrate for **thickness, density, and interface roughness**,
seeded by the existing fringe analysis, with parameter uncertainties and a
model-vs-data overlay. Complements (does not replace) `analyzeFringes`.

## Problem

`xrdc.xrr.analyzeFringes` extracts a thickness from Kiessig-fringe *spacing* only.
That method assumes a single smooth uniform layer with well-resolved fringes, and
returns no density and no roughness; it degrades for thin layers (too few fringes),
rough films (damped fringes), and gives thickness alone. The reference method —
used by GenX and REFL1D — is a slab-model fit of the full reflectivity curve via
the **Parratt recursion** with **Névot–Croce** interface roughness, which recovers
thickness, density (electron-density / SLD), and roughness together, with proper
error bars. This feature adds that, validated on the lab's existing XRR `.txt`
files (S08/S11/S25 etc., all present locally).

This is feature B of a two-feature plan (B: XRR fitting — this spec; A: RSM
strain/composition — separate later spec). See `docs/research/2026-06-12-hrxrd-feature-research.md` Area 2.

## Scope (v1)

- **Single film on a substrate.** One finite layer + a semi-infinite substrate,
  with roughness on the top (film/ambient) and buried (film/substrate) interfaces.
  The Parratt engine is written for a general N-layer stack; v1 only *exposes* the
  single-film case (multilayer/superlattice is a later extension, not a rewrite).
- **Cu Kα only** (λ from the scan; the embedded scattering table is at Cu Kα energy).
- **Named-material + fitted density** optical constants.
- Local optimizer (Levenberg–Marquardt) seeded from the fringe analysis, with a
  derivative-free fallback. No global optimizer in v1.

## Components

### 1. Optical constants — `xrdc.xrr.opticalConstants` (new)

```matlab
[delta, beta] = xrdc.xrr.opticalConstants(material, density, lambda)
```

Computes the dispersion δ and absorption β of the complex refractive index
`n = 1 − δ + iβ` from a material's chemical formula and mass density:

- δ = (rₑ λ² / 2π) · Σᵢ Nᵢ f₁ᵢ
- β = (rₑ λ² / 2π) · Σᵢ Nᵢ f₂ᵢ

where rₑ = 2.8179403×10⁻⁵ Å (classical electron radius), Nᵢ is the number density
of element i (atoms/Å³) computed from the mass density ρ (g/cm³), formula
stoichiometry, and molar mass via Nᵢ = (ρ·N_A/M)·xᵢ with M = Σ xᵢ Aᵢ, and (f₁ᵢ, f₂ᵢ)
are the atomic scattering factors at the X-ray energy (f₁ ≈ Z + f′, f₂ = f″).

- `material`: a name resolvable in the materials database (carries `formula`), or a
  formula string directly (e.g. `"SrTiO3"`).
- `density`: mass density in g/cm³.
- `lambda`: wavelength in Å (used to confirm the table energy; v1 asserts Cu Kα
  within tolerance and errors otherwise — `xrdc:xrr:unsupportedEnergy`).

**New data file:** `+xrdc/+data/atomicScattering.json` — per-element (f₁, f₂) at Cu
Kα (8.0478 keV) for the elements the group uses (at least O, Al, Ti, Sr, Zr, Ru, La,
W, Pb, Pt; extensible). Source values from the Henke/CXRO tables (cited in the file
and in SCIENTIFIC_ASSUMPTIONS). Molar masses Aᵢ from a small embedded element table
(or the same file).

**Materials DB extension:** `+xrdc/+data/materials.json` gains `formula` (string)
and `densityBulk` (g/cm³) fields on each entry, and new entries for substrate/film
materials needed by XRR but absent today: **LaAlO₃ (LAO), TiO₂ (rutile), WO₃**.
Existing consumers (`loadMaterials`, `identifyMaterial`) must continue to work —
new fields are additive; `loadMaterials` returns them but `identifyMaterial` ignores
them.

### 2. Forward model — `xrdc.xrr.reflectivityModel` (new)

```matlab
R = xrdc.xrr.reflectivityModel(twoTheta, layers, lambda, options)
```

Specular reflectivity R(θ) for a stratified medium via the **Parratt recursion**.
`twoTheta` is the scan's 2θ (degrees); the grazing angle is θ = 2θ/2.

- `layers`: struct array, top→bottom, each with `.material`, `.density` (g/cm³),
  `.thickness` (nm; ignored for the substrate), `.roughness` (nm, the roughness of
  that layer's *top* interface). The last element is the semi-infinite substrate.
  Ambient above is vacuum (δ=β=0).
- Per layer, n_j = 1 − δ_j + iβ_j via `opticalConstants`. Vertical wavevector
  k_{z,j} = (2π/λ)·√(n_j² − cos²θ). Interface Fresnel coefficient
  r_{j,j+1} = (k_{z,j} − k_{z,j+1})/(k_{z,j} + k_{z,j+1}), damped by Névot–Croce
  roughness: r′ = r·exp(−2 k_{z,j} k_{z,j+1} σ_j²). Recursion from the substrate
  (X_substrate = 0) upward:
  X_j = (r′_{j,j+1} + X_{j+1} e^{2 i k_{z,j+1} d_{j+1}}) / (1 + r′_{j,j+1} X_{j+1} e^{2 i k_{z,j+1} d_{j+1}}).
  Reflectivity R = |X_top|².
- `options.Footprint` (default `true`): low-angle beam-footprint/illumination
  correction R ← R · min(1, sin θ / sin θ_fp), where θ_fp is the sample-fill angle
  (a model parameter, see fitter). Off → no correction.

Pure function (no instrumental scale/background here — those live in the fitter's
objective so the forward model stays a clean physics kernel).

### 3. Fitter — `xrdc.xrr.fitReflectivity` (new)

```matlab
result = xrdc.xrr.fitReflectivity(scan, options)
```

Fits the single-film model to an XRR scan. Pipeline:

1. **Seed from the existing analysis.** Run `analyzeFringes` for the thickness seed
   and `findCriticalEdge` for θc → δ → density seed (invert the δ formula at the
   film's nominal bulk composition). Roughness seeds small (default 0.3 nm); scale
   from the post-edge plateau; background from the high-angle tail. The film and
   substrate materials are supplied by the caller (`options.Film`,
   `options.Substrate`), defaulting to a reasonable pair if omitted.
2. **Objective in log space** (XRR spans ~6 decades):
   minimize Σ_θ [log10(max(I_obs,1)) − log10(scale·R_model(θ;p) + background)]²
   over the analysis window [θc+buffer, UpperBound] (reusing `analyzeFringes`
   windowing conventions).
3. **Free parameters (single-film):** film thickness, film density, film roughness,
   substrate roughness, scale, background, and (if `Footprint=true`) the footprint
   angle θ_fp. Bounds: thickness>0, 0<density<1.3·bulk, roughness≥0, scale>0, bg≥0.
4. **Optimizer:** `lsqnonlin` (Levenberg–Marquardt / trust-region-reflective with
   bounds) when the Optimization Toolbox is present; **`fminsearch` (Nelder–Mead)**
   fallback otherwise, with positivity enforced by parameter reparameterization
   (fit log of positive params) and a penalty for out-of-bounds. Detection mirrors
   the existing `fitPeak` toolbox-fallback pattern.
5. **Uncertainty:** parameter covariance from the Jacobian (C = σ²(JᵀJ)⁻¹) → ±σ per
   parameter, reusing the approach already in `fitPeak`. Report `NaN` SE (with a
   flag) when the Jacobian is rank-deficient.

**Output `result` (struct):** `.thicknessNm ± .thicknessSeNm`, `.densityGcc ±
.densitySeGcc` (and `.densityFraction` = fitted/bulk), `.filmRoughnessNm`,
`.substrateRoughnessNm` (each ± SE), `.scale`, `.background`, `.footprintDeg`,
`.rSquared`/`.chiSq`, `.modelCurve` (R on the scan's 2θ grid), `.residuals`,
`.window` ([θ_lo θ_hi]), `.converged` (logical), `.method` ("lm"|"neldermead"),
`.seed` (the initial guesses). Error IDs `xrdc:xrr:*`.

### 4. GUI — `xrdcApp.m` XRR panel

Add a **"Fit model"** action to the `xrr` parameter panel (alongside the existing
fringe controls): film/substrate material dropdowns (populated from `materials.json`
entries with `formula`+`densityBulk`), a button that runs `fitReflectivity` seeded by
the current fringe analysis, overlays `result.modelCurve` on the existing log XRR
plot, and writes the fitted parameters ± errors to the Results area. Follows the
established `runXRR` / `onIdentifyMaterial` button pattern and dark-theme conventions.
The existing fringe-spacing readout stays; the fit is additive.

### 5. Demo + docs

- `examples/demoXrrFit.m` — load an XRR scan, fit, print parameters, plot overlay.
- `docs/USER_GUIDE.md` — XRR section: fringe-spacing vs full-fit, when to use which,
  how to call `fitReflectivity` and read the result.
- `docs/FEATURES.md` — status line under the XRR section.
- `docs/SCIENTIFIC_ASSUMPTIONS.md` — new entry: Parratt recursion + Névot–Croce
  assumptions; the atomic-scattering table source and its fixed Cu Kα energy;
  the footprint-correction model; log-space objective; Jacobian-error caveats and
  the XRR phase-problem / fit non-uniqueness note (recommend sanity-checking against
  the fringe-spacing thickness).

## Data flow

scan → `analyzeFringes` + `findCriticalEdge` (seeds) → build initial single-film
layer model (caller-supplied film/substrate materials) → `fitReflectivity`
(Parratt forward via `reflectivityModel`, log-space objective, LM or Nelder–Mead) →
params ± σ + model overlay → GUI/results.

## Error handling

- Unknown material / missing `formula` or `densityBulk` → `xrdc:xrr:unknownMaterial`
  (listing valid names) or `xrdc:xrr:missingMaterialData`.
- Non-Cu-Kα wavelength → `xrdc:xrr:unsupportedEnergy` (v1 table is Cu Kα only).
- No critical edge found → warn (`xrdc:xrr:noEdge`), fall back to a default density
  seed and proceed.
- Optimizer non-convergence → return best-so-far with `.converged = false`
  (not an error).
- Optimization Toolbox absent → Nelder–Mead fallback with a one-time note in the
  result `.method`.
- Scan too short / missing fields → existing `xrdc:xrr:*` guards.

## Testing

`tests/testXrr.m` additions (functiontests, per repo conventions):

- **Optical constants known-answer:** δ/β for SrTiO₃ at Cu Kα → grazing critical
  angle θc = √(2δ) ≈ 0.28° (i.e. 2θc ≈ 0.55°, matching the ~0.55° edge `findCriticalEdge`
  reports on the S-series XRR). Vacuum → δ=β=0. (Convention pinned: the scan's
  `twoTheta` is 2θ; all internal angles use grazing θ = 2θ/2.)
- **Forward-model limits:** bare substrate (no film) reproduces the Fresnel
  reflectivity and a smooth critical edge; a single film produces Kiessig fringes
  whose 2θ spacing matches the analytic Kiessig relation for the set thickness (the
  same relation `thicknessFromFringes` inverts); increasing roughness monotonically
  damps fringe contrast (Névot–Croce).
- **Round-trip recovery:** simulate a single film {d, ρ, σ} via `reflectivityModel`,
  add Poisson noise (fixed `rng`), fit → recover thickness within ~5%, density and
  roughness within tolerance; SEs finite and positive.
- **Seeding:** the fitter's thickness seed equals `analyzeFringes` on a synthetic
  scan; fit improves (lower residual) over the seed.
- **Fallback path:** force the Nelder–Mead branch → round-trip still recovers params.
- **Error paths:** unknown material, non-Cu-Kα λ.
- **Real data (gated `assumeTrue(isfile(...))`, files present locally so they RUN):**
  fit S08 / S11 / S25 XRR `.txt` → thickness within ~10% of the fringe-spacing value
  (S25 ≈ 40 nm); fitted density within a plausible band of bulk; roughness ≥ 0 and
  physically small; `.converged = true`.

## Out of scope (v1)

- Multilayer / superlattice fitting (engine is N-layer; UI exposes single film).
- MCMC / Bayesian posterior uncertainty (Jacobian SEs only).
- Non-Cu-Kα energies; anomalous-dispersion fine structure beyond the tabulated Cu Kα
  f′/f″.
- Graded / SLD-profile (sub-slab) layers; resolution smearing.
- Global optimizer (differential evolution) — local LM from the fringe seed only.
