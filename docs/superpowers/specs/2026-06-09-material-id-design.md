# Design: Material identification from (00l) peak positions

**Date:** 2026-06-09
**Status:** Approved (Shawn, via design review in session)
**Feature:** From a θ-2θ scan of an epitaxial film, identify which material(s) produced
the observed peaks (PZT, SrTiO3, PbTiO3, SrRuO3) and report derived properties.

## Problem

The Paik group grows epitaxial perovskite films (PZT, PTO on STO; SRO electrodes) and
wants the toolkit to answer "what material is this?" from the (00l) peak positions of a
2θ scan — the standard material-ID feature of XRDC-class software. Peaks shift with
epitaxial strain, multiple phases coexist in one scan (substrate + film + electrode),
and PTO vs PZT bulk c values are nearly identical, so naive nearest-reference matching
gives confident wrong answers.

## Approach (chosen: harmonic-series clustering + strain-aware matching)

Rejected alternatives:

- **Per-peak nearest-reference lookup** — a single strained peak cannot distinguish
  PTO(001) from PZT(001); per-peak voting without series consistency is fragile.
- **Full-pattern simulation + cross-correlation** — needs intensity models we don't
  have; overkill for 4 materials in (00l) geometry.

Key physics constraint that shaped the design: **bulk-c matching fails for coherently
strained films.** PTO bulk c (4.152 Å) and MPB-region PZT bulk c (~4.146 Å) differ by
0.006 Å — unresolvable. But pseudomorphic on STO, PZT 52/48 (a_bulk ≈ 4.036 Å) sits
under ~−3.2% in-plane compression, expanding c to ≈ 4.26 Å, while PTO is nearly
lattice-matched and stays ≈ 4.15 Å. Matching must therefore compare the measured c
against the **predicted strained c**, not bulk c. Real films sit anywhere between fully
strained and fully relaxed, so candidates are scored against the interval
[c_bulk, c_pred], not a point.

## Components

### 1. Database — `+xrdc/+data/materials.json` (new)

`substrates.json` (plot overlay) is untouched. Schema per material:

```json
{
  "name": "PbTiO3",
  "aliases": ["PTO"],
  "system": "tetragonal",
  "a": 3.904, "c": 4.152,
  "role": "film",                  // "substrate" | "film" | "both"
  "elastic": { "c13": 0.0, "c33": 0.0 },     // GPa, illustrative — real values cited via refs; cubic entries use {"nu": ...}
  "refs": "literature citation per value"
}
```

Entries: **SrTiO3** (role `both`, cubic, a = 3.905, ν ≈ 0.232), **PbTiO3**,
**PZT** (tetragonal-side composition model: piecewise-linear a(x), c(x) vs Zr fraction
x, valid x ≲ 0.52, with refs), **SrRuO3** (pseudocubic a ≈ 3.93; included because
validated real data — Tushar S25, SRO on STO — gives a free known-answer test).
All numeric values get literature citations in `refs` (confirm exact constants during
implementation; the schema, not the constants, is the design commitment).

Build: add `materials.json` to the bundled-resource list in `build/buildStandalone.m`
(currently `["substrates.json", "xrayLines.json"]`).

### 2. Core function — `xrdc.lattice.identifyMaterial`

```matlab
R = xrdc.lattice.identifyMaterial(twoTheta, lambda, Substrate="SrTiO3", ...)
```

Accepts a vector of peak 2θ positions (degrees) — optionally with intensities/FWHM for
artifact filtering and series quality — or a scan struct (runs `xrdc.peaks.findPeaks`
internally). Pipeline:

1. **Artifact filter.** For each strong peak, compute where its Cu Kβ (1.3922 Å) and
   W Lα (1.4763 Å) ghost reflections would land (same d, different λ); remove weak
   peaks within tolerance of those positions before grouping. Skipped when no
   intensities are provided (with a warning, since strong/weak cannot be ranked).
2. **Substrate confirmation — input, not discovery.** The user declares the substrate
   (default `"SrTiO3"`). Predict its (00l) series from the database, claim observed
   peaks within 0.1%, and assign them. If the series is absent, warn
   (`xrdc:lattice:substrateNotFound`) and continue treating all peaks as unassigned.
3. **Harmonic grouping.** Remaining peaks are clustered into consistent-c series: each
   peak is hypothesised as order l = 1…4 of some c; peaks agreeing within tolerance
   join the series. **c-vs-c/2 guard:** when only even orders are matched (e.g. 002 +
   004 with a weak/absent 001), evaluate both the c and c/2 hypotheses; prefer by fit
   tightness and database plausibility, and flag the ambiguity in the output rather
   than silently choosing.
4. **Per-series refinement.** Least-squares fit of c over all matched orders (higher
   orders dominate precision), with uncertainty from the fit residuals.
5. **Strain-aware naming.** Per database candidate:
   - ε∥ = (a_sub − a_bulk) / a_bulk  (a_sub from the declared substrate entry)
   - ε⊥ = −[2ν/(1−ν)]·ε∥  (cubic)  or  −2(c₁₃/c₃₃)·ε∥  (tetragonal)
   - c_pred = c_bulk·(1 + ε⊥)
   - Score measured c against the interval [c_bulk, c_pred] (tolerance-padded).
   - **Always return the ranked candidate set** within the window — never a lone
     winner. Flag when the runner-up score is close (ambiguous ID).
   - Relaxation fraction r = (c_meas − c_pred) / (c_bulk − c_pred), clamped reporting
     outside [0,1] as "beyond model".
6. **PZT composition.** When PZT is a candidate, estimate Zr fraction x by inverting
   the strain-corrected c(x) model. Output always carries the caveat that composition
   from (00l) alone is strain-confounded and needs in-plane a (RSM or asymmetric
   reflection) to deconvolve.

**Output `R`:** struct with `.series` — a table, one row per identified series:
candidates (ranked sub-table: name, score, strain vs bulk, strain vs pseudomorphic,
relaxation fraction), c_meas ± σ, matched orders/labels, role (substrate/film),
flags (ambiguous, c/2-ambiguous, ghost-filtered count), and for PZT the x estimate +
caveat string; plus `.unassigned` (peaks no series claimed) and `.substrate`
(confirmation result). Error IDs `xrdc:lattice:*`. No new toolbox dependencies.

### 3. GUI integration — `xrdcApp.m`

"Identify Material…" action available for `twoThetaOmega` scans: substrate dropdown
(default SrTiO3, populated from `role` ∈ {substrate, both} entries), runs the pipeline
on the current detected peaks, annotates matched peaks on the plot with material +
(00l) labels, and shows the ranked results in the Analysis panel's Results area
(the app-wide pattern; nested ranked candidates don't fit a uitable), including the PZT
composition caveat line verbatim. Follows the existing dark-theme/white-plot-card
conventions.

## Error handling

- No peaks supplied / none detected → `xrdc:lattice:noPeaks`.
- Unknown substrate name → `xrdc:lattice:unknownMaterial` listing valid names.
- Substrate series not found in scan → warning, not error (scan may be film-only
  range); proceed with naming but mark substrate unconfirmed.
- Series with no database candidate in window → reported as "unidentified" with its
  refined c, never dropped silently.

## Testing

`tests/testIdentify.m` (functiontests), per repo conventions:

- **Synthetic known-answer:** STO substrate + pseudomorphic PTO film (peaks generated
  from c_pred) → PTO ranked first, PZT distinguishable; relaxed PZT 52/48 → PZT first
  with r ≈ 1 behaviour at the bulk end of the interval.
- **Ghost filter:** inject Cu Kβ ghost of the STO 002 → it must not seed a series.
- **c/2 guard:** series with only 002 + 004 present → correct c recovered, ambiguity
  flag set.
- **Ambiguity:** measured c equidistant between two candidates → both reported, flag set.
- **Error paths:** no peaks, unknown substrate.
- **Real data:** S25 SRO/STO scan → SRO identified for the film series; gated with
  `assumeTrue(isfile(...))`.

## Out of scope (v1)

- In-plane parameter / RSM-based strain–composition deconvolution (the caveat points
  users there).
- Non-(00l) orientations and polycrystalline ID.
- Intensity-based scoring (structure factors).
