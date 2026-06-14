# RSM Strain & Composition (Feature A) — Design

**Date:** 2026-06-13
**Status:** Approved (brainstorming complete; ready for implementation plan)
**Author:** Shawn + Claude (Opus 4.8)

## Goal

From an asymmetric reciprocal-space map (RSM), recover a film's in-plane
lattice parameter `a∥` and out-of-plane parameter `a⊥` (= c), decompose the
biaxial strain, recover the *relaxed* (strain-free) pseudocubic parameter
`a₀`, report the degree of relaxation, and — for alloy films with a Vegard
anchor table (PZT) — map `a₀ → composition`. This closes the (00l)-only
degeneracy: a symmetric scan gives only `d⊥`, which entangles composition
and epitaxial strain (the PTO-vs-PZT mis-ranking pinned by
`testIdentifyS11KnownMisrankPtoVsPzt` / `...S31...`). One asymmetric RSM with
`l≠0` and `h²+k²≠0` yields both `a∥` and `c` from a single film peak.

## Scope decisions (from brainstorming)

1. **Peak input — auto-detect in the 2D map.** One call ingests the RSM,
   auto-finds the substrate peak (global intensity max) and the film peak
   (masked secondary max), refines by centroid. A material-independent
   physics core underneath takes peak coordinates directly (for testing).
2. **Reflection — caller specifies `Reflection=[h k l]`.** Substrate and film
   share it (standard for epitaxial RSM). Explicit, zero ambiguity; the
   module independently predicts the substrate peak from the declared
   substrate lattice + this hkl and flags a large mismatch.
3. **Output — strain + lattice numbers for a *declared* film.** No candidate
   ranking, no `identifyMaterial` wire-up this iteration. The module returns
   `a∥, a⊥, a₀, c, ε∥, ε⊥, relaxation, x(PZT)`; the user reads off `(a₀, c)`
   and judges PTO vs PZT. Ranking/ID integration is a clean follow-on once
   real asymmetric perovskite data exists to validate the disambiguation.
4. **Real-data validation — TiO₂ substrate known-answer + add PtO₂.** The
   real PtO₂/TiO₂ 112 RSM validates the geometry core by recovering TiO₂'s
   known a=4.593, c=2.959 from the substrate peak (proven during
   brainstorming: detected at kPar 0.3078 / kPerp 0.6758 → a 4.5946 /
   c 2.9596). A caveated PtO₂ entry is added so `Film="PtO2"` runs the film
   path end-to-end on the real file. Strain decomposition and PZT Vegard are
   additionally proven by synthetic injected-truth known-answer tests.

## Conventions

The existing `xrdc.rsm.toReciprocalSpace` returns `(kPar, kPerp)` in the
**1/d convention**: `|k| = 2 sinθ / λ = 1/d_hkl` (units Å⁻¹), where for a
symmetric (00l) reflection `kPar=0` and `kPerp = l/c`. This module keeps that
convention throughout (no 2π factor). Inversion for a tetragonal/cubic film
at reflection `[h k l]`:

```
a∥ = √(h² + k²) / kPar        (in-plane lattice parameter)
a⊥ = l / kPerp                (out-of-plane parameter, = c)
```

Validated on the real TiO₂ 112 substrate peak to 4 significant figures.

## Physics

**Biaxial strain decomposition** (standard cubic-reference relations; research
report §AREA 1, Birkholz Ch. 7, arXiv:1305.0714):

```
a₀  = ((1−ν)·a⊥ + 2ν·a∥) / (1+ν)
ε⊥  =  (1−ν)/(1+ν) · (a∥ − a⊥) / a₀
ε∥  = −2ν/(1+ν)   · (a∥ − a⊥) / a₀
```

Implemented in terms of the **elastic factor** `f` rather than ν directly, so
both elastic models in `materials.json` work:
- `elastic: {nu}`        → `f = 2ν/(1−ν)`
- `elastic: {c13, c33}`  → `f = 2·c13/c33`

Mirrors the factor `identifyMaterial` already uses for its pseudomorphic
prediction; the shared computation is extracted to one helper.

**Degree of relaxation** (0 = pseudomorphic / fully strained to substrate,
1 = fully relaxed to its own a₀):

```
relaxation = (a∥_film − a∥_substrate) / (a₀ − a∥_substrate)
```

where `a∥_substrate` is the **measured** substrate in-plane parameter from the
same RSM (`√(h²+k²)/kPar` of the substrate peak), not a tabulated bulk value —
the substrate is the unstrained internal reference, so this also absorbs any
shared zero-offset. `pseudomorphic` is true when `relaxation` is within
tolerance of 0.

**Composition (PZT only)** — invert the tetragonal-side Vegard anchor table
`composition.{x,a}` in `materials.json` to map the *relaxed* `a₀ → x`. Uses
relaxed `a₀`, never the clamped `c`. Returns `[]` for films without a
`composition` block. Documented caveat: Vegard is piecewise-valid within one
phase field and breaks at the MPB.

## Architecture — files

### Create

- **`+xrdc/+rsm/biaxialStrain.m`** — pure function
  `[a0, epsPar, epsPerp] = biaxialStrain(aPar, aPerp, factor)`. No I/O, no
  material lookup. The closed-form relations above.

- **`+xrdc/+lattice/elasticFactor.m`** — `f = elasticFactor(material)` where
  `material` is a `loadMaterials` struct. `{nu}` → `2ν/(1−ν)`;
  `{c13,c33}` → `2·c13/c33`. Errors `xrdc:lattice:noElastic` if neither.

- **`+xrdc/+rsm/findRsmPeaks.m`** — material-independent auto-finder.
  `peaks = findRsmPeaks(kPar, kPerp, intensity, Name=Value)`. Returns
  `peaks.substrate{kPar,kPerp,intensity,found}` (global max + intensity-
  weighted centroid in a small window) and `peaks.film{...,found}` (strongest
  local max outside a mask of radius `MaskRadius` around the substrate, above
  `NoiseFactor`×background; `found=false` + no error if none). Options:
  `MaskRadius` (default in Å⁻¹), `NoiseFactor`, `CentroidWindow`.

- **`+xrdc/+rsm/analyzeStrainRSM.m`** — one-call entry point.
  `R = analyzeStrainRSM(rsm, Substrate=, Film=, Reflection=, Name=Value)`.
  `rsm` = folder/file path (→ `loadAreaScan`) or pre-loaded slice struct
  array. Pipeline: per-slice `toReciprocalSpace` → assemble cloud →
  `findRsmPeaks` → hkl inversion → `biaxialStrain` (film) → `a₀`,
  relaxation → Vegard (if PZT) → predict + check substrate peak. Required:
  `Reflection` (1×3). Defaults: `Substrate="SrTiO3"`. `Film` required (must be
  in `materials.json`).

### Modify

- **`+xrdc/+data/materials.json`** — add `PtO2` (film). Literature lattice +
  ν, flagged approximate in `refs`. Keep the identical top-level field set
  (struct-array decode constraint).
- **`+xrdc/+lattice/identifyMaterial.m`** — replace its inline elastic-factor
  computation with a call to the new `elasticFactor` helper (behaviour-
  preserving DRY; existing tests must stay green).

### Tests — `tests/testRsm.m` (extend)

- `biaxialStrain`: zero-strain identity (a∥=a⊥ → a₀=a∥, ε=0); a hand-computed
  strained case; ε∥/ε⊥ sign and the −f/1 ratio.
- `elasticFactor`: ν path, c13/c33 path, error path. (May live in
  `testIdentify.m` next to its consumer — implementer's call.)
- `findRsmPeaks`: synthetic two-Gaussian cloud → both centroids within tol;
  single-peak cloud → `film.found=false`; mask excludes the substrate.
- `analyzeStrainRSM` synthetic: assemble slices with injected substrate+film
  at known (a∥,a⊥); recover a∥, a⊥, a₀, relaxation; PZT case recovers x.
- `analyzeStrainRSM` real-data (gated on `isfile`): the 112 PtO₂/TiO₂ RSM →
  substrate a=4.593±tol, c=2.959±tol; film path returns finite a∥/a⊥/
  relaxation with `Film="PtO2"`.

### Docs / examples

- `examples/demoStrainRSM.m` — run on the real 112 RSM, print the report.
- `docs/USER_GUIDE.md` — new workflow subsection.
- `docs/SCIENTIFIC_ASSUMPTIONS.md` — new section: 1/d convention, biaxial
  relations, pseudocubic approximation caveat, relaxation definition,
  approximate PtO₂ entry, Vegard/MPB caveat.
- `docs/FEATURES.md` — mark RSM strain/composition done.

## Error handling

- Missing `Reflection` → `arguments`-block / explicit error.
- Unknown `Film`/`Substrate` → propagate `xrdc:lattice:unknownMaterial`.
- Film peak not found → `R.film.found=false`, substrate-only result, flag set,
  no error (mirrors the "never silently fail" pattern: surfaced, not hidden).
- Substrate peak far from predicted → `flags` includes `substrateOffPrediction`
  with the percent offset (does not abort; the caller may have a real offset).
- `l=0` or `h=k=0` reflection (no out-of-plane / no in-plane info) →
  `xrdc:rsm:badReflection` (an asymmetric reflection is required).

## Test plan summary

Every new numerical function gets a known-answer test (CLAUDE.md rule). The
geometry core is validated against real TiO₂ truth; strain decomposition and
composition against synthetic injected truth. Full `runtests` green before
each commit; gated real-data test skips cleanly when the file is absent.

## Out of scope (future)

- Candidate ranking by (a, c) and `identifyMaterial` RSM integration (the
  *automatic* PTO-vs-PZT fix) — needs real asymmetric perovskite data first.
- Combining symmetric 002 + asymmetric reflection for higher precision.
- Orthorhombic/monoclinic distortion via the φ = 0/90/180/270° family.
- Interactive / click-to-pick peak selection.
- Refraction correction at grazing incidence/exit on the asymmetric peak.
