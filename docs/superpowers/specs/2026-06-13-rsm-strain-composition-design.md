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

**Biaxial strain decomposition** (standard cubic-reference relations, derived
from Hooke's law with σ_zz = 0 at the free surface; research report §AREA 1,
Birkholz Ch. 7, arXiv:1305.0714):

```
a₀  = ((1−ν)·a⊥ + 2ν·a∥) / (1+ν)   = (a⊥ + f·a∥) / (1+f)
ε∥  = +(1−ν)/(1+ν) · (a∥ − a⊥) / a₀     ← in-plane
ε⊥  = −2ν/(1+ν)   · (a∥ − a⊥) / a₀     ← out-of-plane
```

**Sign/label discipline (this block has been corrected once — do not "tidy"
it).** Physical check: a film under in-plane *compression* has a∥ < a₀ and, by
Poisson, a⊥ > a₀, so (a∥ − a⊥) < 0 → ε∥ < 0 (compressed in-plane ✓),
ε⊥ > 0 (expanded out-of-plane ✓). The invariant the test anchors to is
`ε⊥/ε∥ = −2ν/(1−ν) = −f` — independent of the (a∥−a⊥) magnitude, so it
survives any common-source numeric error. `a₀` is computed *independently* of
ε (so relaxation and PZT composition, which depend only on a₀, are correct
regardless of the ε labels — but the reported ε values go straight into a
methods section, so the labels must be right).

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
tolerance of 0. Note `a₀` (hence relaxation) depends on ν, so an inaccurate ν
biases relaxation *and* both strain components — see the ν-sensitivity caveat.

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
  local max outside a mask around the substrate, above `NoiseFactor`×
  background; `found=false` + no error if none). Options: `MaskRadius`,
  `NoiseFactor`, `CentroidWindow`.
  **Hardening (the auto-finder is weakest exactly where the science is
  hardest).** The feature exists to resolve the *near-degenerate* case — film
  peak close to the substrate — which is precisely when a fixed mask can
  swallow the film peak. Two guards: (1) `MaskRadius` defaults relative to the
  substrate–film separation scale (data-driven), not a fixed Å⁻¹ constant;
  (2) when the detected film peak lies within `1.5 × MaskRadius` of the
  substrate, set a `filmNearSubstrate` flag so the caller knows the hard case
  is in play and the result needs a human glance. Also flag `filmNotBrighter`
  if no clear secondary maximum is found above background — substrate-as-
  brightest fails for a thick/strongly-diffracting film, and silently
  returning the substrate twice would be the worst outcome.

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

- `biaxialStrain`: anchor to assertions the label-swap **cannot survive**, not
  to numbers derived from the formula block (that would let a swapped impl and
  a swapped expected-value pass together — the spec-and-test-from-one-source
  trap):
  - zero-strain identity (a∥=a⊥ → a₀=a∥, ε∥=ε⊥=0);
  - **physical sign** — inject a known in-plane-*compressive* case (a∥ < a⊥)
    and assert `ε∥ < 0` AND `ε⊥ > 0` (the swap flips both);
  - **reconstruction identity** — assert `a∥ ≈ a₀·(1+ε∥)` and
    `a⊥ ≈ a₀·(1+ε⊥)` to round-off (fails immediately if labels are swapped,
    regardless of where the expected numbers came from);
  - **invariant** — assert `ε⊥/ε∥ ≈ −f` for a nonzero-strain case.
- `elasticFactor`: ν path, c13/c33 path, error path. (May live in
  `testIdentify.m` next to its consumer — implementer's call.)
- `findRsmPeaks`: synthetic two-Gaussian cloud → both centroids within tol;
  single-peak cloud → `film.found=false`; mask excludes the substrate;
  **close-pair cloud** (film within 1.5×MaskRadius of substrate) → still
  resolves both AND sets `filmNearSubstrate`; **substrate-only / no clear
  secondary** → `filmNotBrighter` flag rather than returning the substrate
  twice.
- `analyzeStrainRSM` synthetic: assemble slices with injected substrate+film
  at known (a∥,a⊥); recover a∥, a⊥, a₀, relaxation; PZT case recovers x.
- `analyzeStrainRSM` real-data (gated on `isfile`): the 112 PtO₂/TiO₂ RSM →
  substrate a=4.593±tol, c=2.959±tol; film path returns finite a∥/a⊥/
  relaxation with `Film="PtO2"`.

### Docs / examples

- `examples/demoStrainRSM.m` — run on the real 112 RSM, print the report.
- `docs/USER_GUIDE.md` — new workflow subsection.
- `docs/SCIENTIFIC_ASSUMPTIONS.md` — new section covering every item in
  **Caveats & known limitations** below, verbatim in intent. The
  validation-boundary sentence and the pseudocubic/Vegard definitional
  caveat are required, not optional.
- `docs/FEATURES.md` — mark RSM strain/composition done, with the explicit
  note that strain/composition outputs are real-data-unvalidated this
  iteration (geometry core only).

## Error handling

- Missing `Reflection` → `arguments`-block / explicit error.
- Unknown `Film`/`Substrate` → propagate `xrdc:lattice:unknownMaterial`.
- Film peak not found / no clear secondary → `R.film.found=false`, substrate-
  only result, `filmNotBrighter` flag, no error (mirrors the "never silently
  fail" pattern: surfaced, not hidden — and never returns the substrate as if
  it were the film).
- Film peak within `1.5 × MaskRadius` of the substrate → `filmNearSubstrate`
  flag (the near-degenerate hard case; result is returned but wants a human
  glance).
- Substrate peak far from predicted → `flags` includes `substrateOffPrediction`
  with the percent offset (does not abort; the caller may have a real offset).
- `l=0` or `h=k=0` reflection (no out-of-plane / no in-plane info) →
  `xrdc:rsm:badReflection` (an asymmetric reflection is required).

## Test plan summary

Every new numerical function gets a known-answer test (CLAUDE.md rule). The
geometry core is validated against real TiO₂ truth; strain decomposition and
composition against synthetic injected truth. Full `runtests` green before
each commit; gated real-data test skips cleanly when the file is absent.

**Validation boundary (state plainly, do not bury).** The real PtO₂/TiO₂ RSM
validates *only the geometry core* — the hkl inversion against the unstrained
substrate peak. Every part of the headline capability (biaxial decomposition,
a₀ recovery, relaxation, composition, the PTO-vs-PZT disambiguation) is
validated *only against synthetic injected truth*, which tests that the
algebra inverts the forward model — not that the forward model matches a real
strained film. Whether ν for PTO/PtO₂ is right, whether the pseudocubic
approximation holds, and whether real peak shapes centroid where expected are
all **untested against ground truth as of this iteration**.

## Caveats & known limitations

These go into `SCIENTIFIC_ASSUMPTIONS.md`. The first two are required.

1. **Strain & composition are real-data-unvalidated this iteration.** See the
   validation boundary above. A reader should not trust a reported ε or x
   number without an independent check until a real strained film with known
   parameters has been run through.
2. **Pseudocubic approximation — and a Vegard definitional-consistency risk.**
   For PTO/PZT the *relaxed* crystal is itself tetragonal (c/a ≈ 1.06), not
   cubic. The biaxial decomposition assumes a single cubic reference a₀, so the
   recovered a₀ is a pseudocubic strain-model *average* that conflates
   spontaneous tetragonality with epitaxial strain — it is **not** a physical
   relaxed lattice constant. Mapping it through Vegard is therefore doubly
   approximate, and carries a concrete bias risk: the composition table's
   tabulated `a` must be defined the *same way* as a₀ (pseudocubic average vs
   true tetragonal a-axis). A table of bulk a-axis values fed a pseudocubic a₀
   yields a clean-looking but systematically biased x. The implementer must
   confirm the `composition.a` anchors and a₀ use a consistent definition (and
   document which). Sits next to the existing MPB piecewise-validity caveat.
3. **ν sensitivity.** a₀, relaxation, and both strain components all depend on
   ν; an inaccurate ν biases all of them, not just strain. PtO₂'s ν is a
   placeholder — flag the dependence.
4. **Tilt assumed zero.** A tilted film (mosaic/miscut) biases kPar and hence
   a∥ (and through it relaxation/composition). Out of scope this iteration; one
   sentence in the assumptions doc so it isn't mistaken for handled.
5. **Centroid bias from the analyzer/CTR streak.** 2D background and the
   crystal-truncation-rod / analyzer streak through the substrate peak can pull
   the intensity-weighted centroid; the substrate centroid is the more
   exposed of the two. Note as a known small systematic.
6. **`materials.json` identical-field tax.** Every new entry must mirror all
   top-level fields (including empty elastic/composition placeholders) for the
   struct-array decode. Tolerated now; a normalizing loader would remove the
   footgun later (out of scope).

## Out of scope (future)

- Candidate ranking by (a, c) and `identifyMaterial` RSM integration (the
  *automatic* PTO-vs-PZT fix) — needs real asymmetric perovskite data first.
- Combining symmetric 002 + asymmetric reflection for higher precision.
- Orthorhombic/monoclinic distortion via the φ = 0/90/180/270° family.
- Interactive / click-to-pick peak selection.
- Refraction correction at grazing incidence/exit on the asymmetric peak.
