# RSM Strain & Composition (Feature A) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** From an asymmetric RSM, auto-detect substrate + film peaks, invert to in-plane/out-of-plane lattice parameters, decompose biaxial strain to a relaxed a₀, report degree of relaxation, and (for PZT) composition — for a declared film.

**Architecture:** Four focused units under `+xrdc/` — a pure strain function (`biaxialStrain`), a shared elastic-factor helper (`elasticFactor`, extracted from `identifyMaterial`), a material-independent 2D peak auto-finder (`findRsmPeaks`), and a one-call orchestrator (`analyzeStrainRSM`) that wires `loadAreaScan`→`toReciprocalSpace`→finder→inversion→strain→composition. Geometry validated against the real TiO₂ substrate peak; strain/composition against synthetic injected truth.

**Tech Stack:** MATLAB R2022b+, `arguments` blocks, `string` identifiers, `functiontests(localfunctions)`, error IDs `xrdc:<subpkg>:<reason>`. Spec: `docs/superpowers/specs/2026-06-13-rsm-strain-composition-design.md`.

**Base MATLAB only — NO Image Processing Toolbox** (CLAUDE.md: no unlicensed-toolbox deps). The peak finder uses `accumarray`, `discretize`, `conv2`, `sub2ind` — NOT `imgaussfilt`, `ordfilt2`, or `imregionalmax`. The Optimization Toolbox is also not needed here. Every algorithm in this plan was prototyped against real and synthetic data before the plan was written; the implementations below are the validated ones — match them.

**Conventions reused (read before starting):**
- `xrdc.rsm.toReciprocalSpace(scan)` → `[kPar, kPerp]` in the **1/d convention**: `|k| = 2 sinθ/λ`, so for reflection `[h k l]`: `a∥ = √(h²+k²)/kPar`, `a⊥ = l/kPerp`.
- `xrdc.lattice.loadMaterials(name)` → struct with `.a .c .elastic .composition .name`.
- `materials.json` entries share an identical top-level field set (struct-array decode).
- Run a single test file during TDD: `matlab -batch "addpath(pwd); addpath('tests'); run(testRsm)"` (substitute `testIdentify` for Task 1). Run the full suite before each commit: `matlab -batch "r = runtests; assert(all([r.Passed]), 'suite not green')"`. The project's `runtests.m` is a zero-arg wrapper over `tests/`.
- **NEVER push.** Commit after each green task. The owner pushes.

---

### Task 1: `elasticFactor` shared helper (extract from `identifyMaterial`, DRY)

**Files:**
- Create: `xrdc-matlab/+xrdc/+lattice/elasticFactor.m`
- Modify: `xrdc-matlab/+xrdc/+lattice/identifyMaterial.m:282-288` (replace local `strainFactor` body with a call)
- Test: `xrdc-matlab/tests/testIdentify.m` (add 1 test)

Background: `identifyMaterial` has a private `strainFactor(e)` (lines 282-288) computing `2ν/(1−ν)` or `2·c13/c33`. `biaxialStrain` (Task 2) needs the same. Extract it to a public function; keep `identifyMaterial` behaviour identical.

- [ ] **Step 1: Write the failing test** — append to `tests/testIdentify.m` (before the final helper block; it's a localfunction so order is free):

```matlab
function testElasticFactorNuAndC13C33(tc)
    % nu path: f = 2 nu / (1 - nu)
    sto = xrdc.lattice.loadMaterials("SrTiO3");   % nu = 0.232
    tc.verifyEqual(xrdc.lattice.elasticFactor(sto), ...
        2*0.232/(1-0.232), 'AbsTol', 1e-12);
    % c13/c33 path: f = 2 c13/c33
    e = struct('elastic', struct('c13', 70, 'c33', 90));
    tc.verifyEqual(xrdc.lattice.elasticFactor(e), 2*70/90, 'AbsTol', 1e-12);
    % error path: neither present
    bad = struct('elastic', struct());
    tc.verifyError(@() xrdc.lattice.elasticFactor(bad), 'xrdc:lattice:noElastic');
end
```

- [ ] **Step 2: Run it, verify it fails** — `matlab -batch "addpath(pwd); addpath('tests'); run(testIdentify)"`. Expected: `testElasticFactorNuAndC13C33` errors (undefined `elasticFactor`).

- [ ] **Step 3: Create the helper** — `+xrdc/+lattice/elasticFactor.m`:

```matlab
function f = elasticFactor(material)
%ELASTICFACTOR  Biaxial strain factor f from a material's elastic block.
%   f = xrdc.lattice.elasticFactor(material)
%
%   f relates out-of-plane to in-plane strain as eps_perp = -f * eps_par for
%   a biaxially strained (001) film. Two supported elastic models:
%     elastic.nu          -> f = 2*nu/(1-nu)     (isotropic Poisson)
%     elastic.{c13,c33}   -> f = 2*c13/c33       (tetragonal single-crystal)
%
%   Input
%     material : struct with an .elastic field (e.g. from loadMaterials).
%
%   Errors: xrdc:lattice:noElastic if neither model is present.

    arguments
        material (1,1) struct
    end
    if ~isfield(material, 'elastic') || ~isstruct(material.elastic)
        error('xrdc:lattice:noElastic', ...
            'Material has no elastic block; cannot compute strain factor.');
    end
    e = material.elastic;
    if isfield(e, 'nu') && ~isempty(e.nu)
        f = 2 * e.nu / (1 - e.nu);
    elseif isfield(e, 'c13') && isfield(e, 'c33') && ~isempty(e.c13) && ~isempty(e.c33)
        f = 2 * e.c13 / e.c33;
    else
        error('xrdc:lattice:noElastic', ...
            'elastic block lacks both nu and {c13,c33}.');
    end
end
```

- [ ] **Step 4: Refactor `identifyMaterial`** — replace the body of its private `strainFactor` (lines 282-288) so it delegates (keeps the call sites unchanged):

```matlab
function f = strainFactor(e)
    f = xrdc.lattice.elasticFactor(e);
end
```

- [ ] **Step 5: Run the new test + the whole identify suite** — `matlab -batch "addpath(pwd); addpath('tests'); r = run(testIdentify); assert(all([r.Passed]), 'testIdentify not green'); disp('OK')"`. Expected: all pass (the 29 existing + 1 new = 30). The refactor is behaviour-preserving, so the S11/S31 characterization tests must still pass unchanged.

- [ ] **Step 6: Full suite + commit**

```bash
matlab -batch "r = runtests; assert(all([r.Passed]), 'suite not green')"
git add "xrdc-matlab/+xrdc/+lattice/elasticFactor.m" "xrdc-matlab/+xrdc/+lattice/identifyMaterial.m" "xrdc-matlab/tests/testIdentify.m"
git commit -m "refactor(lattice): extract elasticFactor helper (DRY for RSM strain)"
```

---

### Task 2: `biaxialStrain` pure physics function

**Files:**
- Create: `xrdc-matlab/+xrdc/+rsm/biaxialStrain.m`
- Test: `xrdc-matlab/tests/testRsm.m` (add tests)

Physics (spec "Physics" block — corrected, do not re-derive): given in-plane `aPar`, out-of-plane `aPerp`, and elastic factor `f`,
`a0 = (aPerp + f·aPar)/(1+f)`, `epsPar = (aPar−a0)/a0`, `epsPerp = (aPerp−a0)/a0`.
This compact form is algebraically identical to the spec's `±` relations and satisfies `epsPerp/epsPar = −f` and the reconstruction identity by construction.

- [ ] **Step 1: Write the failing tests** — append localfunctions to `tests/testRsm.m`:

```matlab
function testBiaxialStrainZeroWhenCubic(tc)
    % a_par == a_perp -> unstrained: a0 = a_par, both strains zero.
    [a0, ep, ez] = xrdc.rsm.biaxialStrain(3.905, 3.905, 0.857);
    tc.verifyEqual(a0, 3.905, 'AbsTol', 1e-12);
    tc.verifyEqual(ep, 0, 'AbsTol', 1e-12);
    tc.verifyEqual(ez, 0, 'AbsTol', 1e-12);
end

function testBiaxialStrainCompressivePhysicalSign(tc)
    % In-plane compression: a_par < a_perp (Poisson pushes c out).
    % Correct labels => eps_par < 0 (compressed) AND eps_perp > 0 (expanded).
    % A label swap flips BOTH signs, so this fails on the swap regardless of
    % where any expected numbers came from.
    f = 0.857;                         % ~ nu=0.30
    [~, epsPar, epsPerp] = xrdc.rsm.biaxialStrain(3.90, 4.05, f);
    tc.verifyLessThan(epsPar, 0, 'in-plane compression must give eps_par < 0');
    tc.verifyGreaterThan(epsPerp, 0, 'Poisson expansion must give eps_perp > 0');
end

function testBiaxialStrainReconstructionAndInvariant(tc)
    % Reconstruction identity: a = a0*(1+eps). Fails immediately if labels
    % are swapped. Invariant: eps_perp/eps_par = -f, independent of magnitude.
    f = 2*0.232/(1-0.232);
    aPar = 3.88; aPerp = 4.02;
    [a0, epsPar, epsPerp] = xrdc.rsm.biaxialStrain(aPar, aPerp, f);
    tc.verifyEqual(a0*(1+epsPar),  aPar,  'AbsTol', 1e-12, 'reconstruct a_par');
    tc.verifyEqual(a0*(1+epsPerp), aPerp, 'AbsTol', 1e-12, 'reconstruct a_perp');
    tc.verifyEqual(epsPerp/epsPar, -f, 'AbsTol', 1e-12, 'eps_perp/eps_par = -f');
end

function testBiaxialStrainRejectsNonPositive(tc)
    tc.verifyError(@() xrdc.rsm.biaxialStrain(0, 4.0, 0.8), 'MATLAB:validators:mustBePositive');
end
```

- [ ] **Step 2: Run, verify failure** — `matlab -batch "addpath(pwd); addpath('tests'); run(testRsm)"`. Expected: the four `testBiaxialStrain*` error (undefined function).

- [ ] **Step 3: Implement** — `+xrdc/+rsm/biaxialStrain.m`:

```matlab
function [a0, epsPar, epsPerp] = biaxialStrain(aPar, aPerp, factor)
%BIAXIALSTRAIN  Relaxed lattice parameter and biaxial strain of a (001) film.
%   [a0, epsPar, epsPerp] = xrdc.rsm.biaxialStrain(aPar, aPerp, factor)
%
%   Closed-form biaxial-strain decomposition for a (001)-oriented film
%   (Hooke's law, sigma_zz = 0 at the free surface). Given the measured
%   in-plane parameter aPar, out-of-plane parameter aPerp, and the elastic
%   factor f (= 2nu/(1-nu) or 2 c13/c33; see xrdc.lattice.elasticFactor):
%
%       a0      = (aPerp + f*aPar) / (1 + f)        relaxed (pseudocubic) param
%       epsPar  = (aPar  - a0) / a0                 in-plane strain  (<0 = compr.)
%       epsPerp = (aPerp - a0) / a0                 out-of-plane strain
%
%   These satisfy epsPerp/epsPar = -f and the reconstruction identity
%   a = a0*(1+eps) exactly. For an intrinsically tetragonal film (PTO/PZT)
%   a0 is the pseudocubic strain-model average, NOT a physical relaxed cubic
%   constant -- see docs/SCIENTIFIC_ASSUMPTIONS.md.
%
%   Inputs are scalar lattice parameters in Angstrom; factor is dimensionless.

    arguments
        aPar   (1,1) double {mustBePositive}
        aPerp  (1,1) double {mustBePositive}
        factor (1,1) double {mustBeNonnegative}
    end
    a0      = (aPerp + factor*aPar) / (1 + factor);
    epsPar  = (aPar  - a0) / a0;
    epsPerp = (aPerp - a0) / a0;
end
```

- [ ] **Step 4: Run, verify pass** — same command as Step 2. Expected: the four `testBiaxialStrain*` pass.

- [ ] **Step 5: Full suite + commit**

```bash
matlab -batch "r = runtests; assert(all([r.Passed]), 'suite not green')"
git add "xrdc-matlab/+xrdc/+rsm/biaxialStrain.m" "xrdc-matlab/tests/testRsm.m"
git commit -m "feat(rsm): biaxialStrain decomposition (relaxed a0 + strain)"
```

---

### Task 3: Add caveated PtO₂ entry to `materials.json`

**Files:**
- Modify: `xrdc-matlab/+xrdc/+data/materials.json` (add one element after the WO3 entry)
- Test: `xrdc-matlab/tests/testIdentify.m` (add 1 test)

PtO₂ is needed so `Film="PtO2"` runs the film path end-to-end on the real RSM. Values are approximate (β-PtO₂, CaCl₂-type, treated pseudo-tetragonally) and flagged as such; the rigorous substrate-side test does not depend on them. **Keep the identical top-level field set** (`name, aliases, system, a, c, role, elastic, composition, refs, formula, densityBulk`).

- [ ] **Step 1: Write the failing test** — append to `tests/testIdentify.m`:

```matlab
function testLoadMaterialsHasPtO2(tc)
    e = xrdc.lattice.loadMaterials("PtO2");
    tc.verifyEqual(string(e.name), "PtO2");
    tc.verifyEqual(string(e.system), "tetragonal");
    tc.verifyGreaterThan(e.a, 4.0);   % ~4.485
    tc.verifyGreaterThan(e.c, 3.0);   % ~3.137
    tc.verifyTrue(isfield(e.elastic, 'nu') && ~isempty(e.elastic.nu));
    % approximate values must be advertised in refs
    tc.verifyTrue(contains(lower(string(e.refs)), "approximate"));
end
```

- [ ] **Step 2: Run, verify failure** — `matlab -batch "addpath(pwd); addpath('tests'); run(testIdentify)"`. Expected: `testLoadMaterialsHasPtO2` fails (unknown material `PtO2`).

- [ ] **Step 3: Add the entry** — in `+xrdc/+data/materials.json`, insert this element immediately after the closing `}` of the `WO3` entry (add a comma after WO3's `}` and paste before the closing `]`):

```json
    ,{
      "name": "PtO2",
      "aliases": ["platinum dioxide", "beta-PtO2"],
      "system": "tetragonal",
      "a": 4.485,
      "c": 3.137,
      "role": "film",
      "elastic": { "nu": 0.30 },
      "composition": null,
      "refs": "APPROXIMATE: beta-PtO2 (CaCl2-type orthorhombic Pnnm a=4.485 b=4.533 c=3.137 A) treated pseudo-tetragonally (a~=b); nu=0.30 placeholder. Used only so Film=PtO2 runs end-to-end; substrate-side RSM validation does not depend on these values.",
      "formula": "PtO2",
      "densityBulk": 11.8
    }
```

- [ ] **Step 4: Run, verify pass + struct-array decode intact** — `matlab -batch "addpath(pwd); addpath('tests'); r = run(testIdentify); assert(all([r.Passed]),'testIdentify not green'); disp('OK')"`. Expected: all pass. (If `loadMaterials` errors on decode, a field set mismatch was introduced — fix the new entry's fields.)

- [ ] **Step 5: Full suite + commit**

```bash
matlab -batch "r = runtests; assert(all([r.Passed]), 'suite not green')"
git add "xrdc-matlab/+xrdc/+data/materials.json" "xrdc-matlab/tests/testIdentify.m"
git commit -m "data(materials): add caveated PtO2 entry for RSM film path"
```

---

### Task 4: `findRsmPeaks` — material-independent 2D auto-finder

**Files:**
- Create: `xrdc-matlab/+xrdc/+rsm/findRsmPeaks.m`
- Test: `xrdc-matlab/tests/testRsm.m` (add tests + a synthetic-cloud helper)

Takes flat `(kPar, kPerp, intensity)` point clouds (the assembled RSM) and returns substrate + film peaks via grid + regional-maxima (see the Algorithm note below — the naive mask approach does NOT separate close peaks). Substrate = brightest regional max; film = next regional max ≥3 grid cells away, above `NoiseFactor`×background; both centroid-refined. Flags: `filmNearSubstrate` when the film is within `NearThreshold` of the substrate (near-degenerate, returned but wants a glance); `filmNotBrighter` (and `film.found=false`) when no clear secondary exists — the substrate is never returned as the film.

- [ ] **Step 1: Write the failing tests + helper** — append to `tests/testRsm.m`. The helper makes a sharp substrate (single crystals are sharp) and a broader film, so the substrate's shoulder never out-shines a nearby film:

```matlab
function [kp, kz, I] = twoGaussianCloud(cSub, cFilm, ampFilm)
%Helper: dense grid cloud, sharp substrate (sigma 0.0025) + broader film
%   (sigma 0.004), at cSub, cFilm (=[kPar kPerp]). Substrate brightest.
    [KP, KZ] = meshgrid(linspace(0.20, 0.40, 200), linspace(0.60, 0.85, 220));
    wS = 0.0025; wF = 0.004;
    I = 1 ...
      + 1000   *exp(-((KP-cSub(1)).^2  + (KZ-cSub(2)).^2 )/(2*wS^2)) ...
      + ampFilm*exp(-((KP-cFilm(1)).^2 + (KZ-cFilm(2)).^2)/(2*wF^2));
    kp = KP(:); kz = KZ(:); I = I(:);
end

function testFindRsmPeaksTwoWellSeparated(tc)
    cSub = [0.256 0.768]; cFilm = [0.244 0.732];   % sep ~0.038
    [kp, kz, I] = twoGaussianCloud(cSub, cFilm, 300);
    p = xrdc.rsm.findRsmPeaks(kp, kz, I);
    tc.verifyTrue(p.substrate.found);
    tc.verifyTrue(p.film.found);
    tc.verifyEqual([p.substrate.kPar p.substrate.kPerp], cSub, 'AbsTol', 5e-3);
    tc.verifyEqual([p.film.kPar p.film.kPerp], cFilm, 'AbsTol', 5e-3);
    tc.verifyFalse(any(p.flags == "filmNearSubstrate"));
end

function testFindRsmPeaksNearDegenerateFlags(tc)
    % Film close to substrate (sep ~0.016 < NearThreshold 0.03) -> still
    % resolved by the grid maxima, but flagged.
    cSub = [0.256 0.768]; cFilm = [0.256 0.752];
    [kp, kz, I] = twoGaussianCloud(cSub, cFilm, 400);
    p = xrdc.rsm.findRsmPeaks(kp, kz, I);
    tc.verifyTrue(p.film.found);
    tc.verifyEqual([p.film.kPar p.film.kPerp], cFilm, 'AbsTol', 5e-3);
    tc.verifyTrue(any(p.flags == "filmNearSubstrate"));
end

function testFindRsmPeaksSubstrateOnly(tc)
    % Single peak -> film not found, flagged, and NOT a copy of substrate.
    cSub = [0.256 0.768];
    [kp, kz, I] = twoGaussianCloud(cSub, [0.30 0.70], 0);   % film amp 0
    p = xrdc.rsm.findRsmPeaks(kp, kz, I);
    tc.verifyTrue(p.substrate.found);
    tc.verifyFalse(p.film.found);
    tc.verifyTrue(any(p.flags == "filmNotBrighter"));
end
```

- [ ] **Step 2: Run, verify failure** — `matlab -batch "addpath(pwd); addpath('tests'); run(testRsm)"`. Expected: the three `testFindRsmPeaks*` error (undefined function).

- [ ] **Step 3: Implement** — `+xrdc/+rsm/findRsmPeaks.m`:

**Algorithm note (validated by prototype on the real 112 RSM + synthetic close/separated/single pairs):** a fixed exclusion mask around the substrate CANNOT separate a film that sits inside the mask — the exact near-degenerate case the feature exists for. Instead: bin the cloud onto a regular grid (max-intensity per cell), lightly smooth, find **regional maxima** (cells ≥ all 8 neighbours), rank by intensity. Substrate = brightest; film = next regional max ≥ 3 grid cells away (so the substrate's own shoulder is never returned as the film). **Base MATLAB only** — `accumarray`/`discretize`/`conv2`/`sub2ind`. Do NOT use `imgaussfilt`, `ordfilt2`, `imregionalmax` or anything from the Image Processing Toolbox (CLAUDE.md: no unlicensed-toolbox deps).

```matlab
function peaks = findRsmPeaks(kPar, kPerp, intensity, options)
%FINDRSMPEAKS  Locate substrate and film peaks in an RSM point cloud.
%   peaks = xrdc.rsm.findRsmPeaks(kPar, kPerp, intensity, Name=Value)
%
%   Material-independent. Bins the scattered (kPar,kPerp,intensity) cloud onto
%   a regular grid (max per cell), lightly smooths, and takes regional maxima.
%   Substrate = brightest regional max (single crystal -> brightest); film =
%   the next-brightest regional max at least 3 grid cells away, above
%   NoiseFactor*median(intensity). Each is refined by an intensity-weighted
%   centroid. Base MATLAB only (no Image Processing Toolbox).
%
%   Name/Value
%     GridBins      grid resolution per axis (default 200).
%     NoiseFactor   a maximum must exceed NoiseFactor*median(intensity)
%                   (default 5).
%     NearThreshold film flagged "filmNearSubstrate" if within this distance
%                   (k-units) of the substrate (default 0.03).
%     RefineWindow  half-width of the centroid window (k-units, default 0.006).
%
%   Output peaks (struct)
%     .substrate .kPar .kPerp .intensity .found
%     .film      .kPar .kPerp .intensity .found
%     .flags     string array: "filmNearSubstrate", "filmNotBrighter"
%
%   Flags (never silently fail): filmNearSubstrate for the near-degenerate
%   hard case (result returned, wants a human glance); filmNotBrighter when no
%   clear secondary maximum exists (film.found=false -- the substrate is never
%   returned as the film).

    arguments
        kPar      (:,1) double
        kPerp     (:,1) double
        intensity (:,1) double
        options.GridBins      (1,1) double {mustBeInteger, mustBePositive} = 200
        options.NoiseFactor   (1,1) double = 5
        options.NearThreshold (1,1) double = 0.03
        options.RefineWindow  (1,1) double = 0.006
    end

    nb = options.GridBins;
    bg = median(intensity);
    kpEdges = linspace(min(kPar),  max(kPar),  nb+1);
    kzEdges = linspace(min(kPerp), max(kPerp), nb+1);
    ip = discretize(kPar,  kpEdges);
    iz = discretize(kPerp, kzEdges);
    valid = ~isnan(ip) & ~isnan(iz);
    G = accumarray([iz(valid) ip(valid)], intensity(valid), [nb nb], @max, 0);

    % light Gaussian smoothing via conv2 (base MATLAB; NO imgaussfilt)
    [xx, yy] = meshgrid(-2:2, -2:2);
    ker = exp(-(xx.^2 + yy.^2)/2); ker = ker / sum(ker(:));
    Gs = conv2(G, ker, 'same');

    % interior regional maxima: cell >= all 8 neighbours (NO ordfilt2)
    cc = Gs(2:end-1, 2:end-1);
    ge = true(size(cc));
    for dr = -1:1
        for dc = -1:1
            if dr == 0 && dc == 0, continue, end
            ge = ge & cc >= Gs(2+dr:end-1+dr, 2+dc:end-1+dc);
        end
    end
    isMax = false(size(Gs));
    isMax(2:end-1, 2:end-1) = ge & cc > options.NoiseFactor*bg;

    kpC = (kpEdges(1:end-1) + kpEdges(2:end))/2;
    kzC = (kzEdges(1:end-1) + kzEdges(2:end))/2;
    dk  = hypot(kpEdges(2)-kpEdges(1), kzEdges(2)-kzEdges(1));
    flags = string.empty(1,0);

    [ri, ci] = find(isMax);
    if isempty(ri)                              % degenerate: fall back to global max
        [mx, im] = max(intensity);
        peaks.substrate = struct('kPar', kPar(im), 'kPerp', kPerp(im), ...
            'intensity', mx, 'found', true);
        peaks.film = struct('kPar', NaN, 'kPerp', NaN, 'intensity', NaN, 'found', false);
        peaks.flags = "filmNotBrighter";
        return
    end
    vals = Gs(sub2ind(size(Gs), ri, ci));
    [vals, ord] = sort(vals, 'descend'); ri = ri(ord); ci = ci(ord);

    subSeed = [kpC(ci(1)), kzC(ri(1))];
    sub = refineCentroid(kPar, kPerp, intensity, subSeed, options.RefineWindow);
    peaks.substrate = struct('kPar', sub(1), 'kPerp', sub(2), ...
        'intensity', vals(1), 'found', true);

    peaks.film = struct('kPar', NaN, 'kPerp', NaN, 'intensity', NaN, 'found', false);
    for m = 2:numel(ri)
        cand = [kpC(ci(m)), kzC(ri(m))];
        if hypot(cand(1)-subSeed(1), cand(2)-subSeed(2)) > 3*dk
            fc = refineCentroid(kPar, kPerp, intensity, cand, options.RefineWindow);
            peaks.film = struct('kPar', fc(1), 'kPerp', fc(2), ...
                'intensity', vals(m), 'found', true);
            if hypot(fc(1)-sub(1), fc(2)-sub(2)) < options.NearThreshold
                flags(end+1) = "filmNearSubstrate"; %#ok<AGROW>
            end
            break
        end
    end
    if ~peaks.film.found
        flags(end+1) = "filmNotBrighter";
    end
    peaks.flags = flags;
end

function c = refineCentroid(kPar, kPerp, intensity, seed, win)
%REFINECENTROID  Intensity-weighted centroid within +/- win of seed (k-units).
    sel = abs(kPar - seed(1)) <= win & abs(kPerp - seed(2)) <= win;
    w = intensity(sel);
    w = max(w - min(w), 0);               % local-background subtract
    if sum(w) <= 0, c = seed; return, end
    c = [sum(kPar(sel).*w), sum(kPerp(sel).*w)] / sum(w);
end
```

- [ ] **Step 4: Run, verify pass** — same command as Step 2. Expected: the three `testFindRsmPeaks*` pass.

- [ ] **Step 5: Full suite + commit**

```bash
matlab -batch "r = runtests; assert(all([r.Passed]), 'suite not green')"
git add "xrdc-matlab/+xrdc/+rsm/findRsmPeaks.m" "xrdc-matlab/tests/testRsm.m"
git commit -m "feat(rsm): findRsmPeaks 2D auto-finder with near-degenerate guards"
```

---

### Task 5: `analyzeStrainRSM` orchestrator + real-data test

**Files:**
- Create: `xrdc-matlab/+xrdc/+rsm/analyzeStrainRSM.m`
- Test: `xrdc-matlab/tests/testRsm.m` (add synthetic + gated real-data tests + a slice-builder helper)

One-call entry: `R = analyzeStrainRSM(rsm, Substrate=, Film=, Reflection=)`. `rsm` is a folder/file path (→ `loadAreaScan`) or a pre-loaded slice struct array. Assembles the q-cloud via `toReciprocalSpace` per slice, runs `findRsmPeaks`, inverts both peaks to `(a∥, a⊥)`, decomposes the film strain, computes relaxation against the **measured** substrate a∥, and (for PZT) composition. Predicts the substrate peak from its lattice and flags a large offset.

- [ ] **Step 1: Write the failing tests + helper** — append to `tests/testRsm.m`. The helper builds twoThetaOmega slices spanning two target k-points and injects two Gaussians **in k-space**, so the peaks land exactly at the chosen `(kPar,kPerp)` regardless of the angle grid:

```matlab
function scans = makeRsmScans(hkl, kSub, kFilm, ampFilm, lambda) %#ok<INUSL>
%Helper: build twoThetaOmega slices whose assembled q-cloud has substrate+film
%   Gaussians exactly at kSub, kFilm (=[kPar kPerp]). Geometry note: for a
%   fixed secondAxis the asymmetric transform traces a RADIAL line at angle
%   phi = secondAxis - 2theta_ctr/2 (constant along the slice), radius
%   (2/lambda)sin(theta). So to put a slice through a target's phi we set
%   secondAxis = phi_target + 2theta_ctr/2; the target lands at that slice's
%   2theta = 2*asin(lambda|k|/2). (hkl is unused here but kept for caller
%   symmetry / documentation.)
    kAll = [kSub; kFilm];
    kmag = hypot(kAll(:,1), kAll(:,2));
    tt   = 2*asin(lambda*kmag/2)*180/pi;            % 2theta per target (deg)
    phi  = atan2(kAll(:,1), kAll(:,2))*180/pi;       % radial angle per target (deg)
    pad  = 2.0;
    ttRange = [min(tt)-pad, max(tt)+pad];
    ttGrid  = linspace(ttRange(1), ttRange(2), 240).';
    ttCtr   = mean(ttRange);                          % = 2theta_ctr the transform uses
    sa      = phi + ttCtr/2;                           % secondAxis hitting each phi
    omRange = [min(sa)-pad, max(sa)+pad];
    omGrid  = linspace(omRange(1), omRange(2), 100);
    wS = 0.002; wF = 0.005;                            % sharp substrate, broad film
    scans = repmat(xrdc.io.emptyScan(), 1, numel(omGrid));
    for j = 1:numel(omGrid)
        s = xrdc.io.emptyScan();
        s.twoTheta = ttGrid; s.secondAxis = omGrid(j);
        s.secondAxisName = "Omega"; s.scanType = "twoThetaOmega"; s.lambda = lambda;
        [kp, kz] = xrdc.rsm.toReciprocalSpace(s);
        I = 1 + 1000  *exp(-((kp-kSub(1)).^2 +(kz-kSub(2)).^2 )/(2*wS^2)) ...
              + ampFilm*exp(-((kp-kFilm(1)).^2+(kz-kFilm(2)).^2)/(2*wF^2));
        s.counts = I;
        scans(j) = s;
    end
end

function testAnalyzeStrainRSMSyntheticPseudomorphic(tc)
    % SrTiO3 [1 0 3] substrate; film pseudomorphic in-plane (a_par = a_sub),
    % tetragonally strained out-of-plane (a_perp = 4.10).
    lambda = 1.5406; hkl = [1 0 3]; hk = hypot(hkl(1), hkl(2));
    aSub = 3.905; aParF = 3.905; aPerpF = 4.10;
    kSub  = [hk/aSub,  hkl(3)/aSub];
    kFilm = [hk/aParF, hkl(3)/aPerpF];
    scans = makeRsmScans(hkl, kSub, kFilm, 350, lambda);
    R = xrdc.rsm.analyzeStrainRSM(scans, Substrate="SrTiO3", Film="PbTiO3", Reflection=hkl);
    tc.verifyEqual(R.substrate.aMeas, aSub,  'AbsTol', 5e-3);
    tc.verifyEqual(R.aPar,  aParF,  'AbsTol', 5e-3);
    tc.verifyEqual(R.aPerp, aPerpF, 'AbsTol', 5e-3);
    % pseudomorphic: a_par == a_sub -> relaxation ~ 0
    tc.verifyLessThan(abs(R.relaxation), 0.1);
    tc.verifyTrue(R.pseudomorphic);
    % a0 matches the pure-function result
    f = xrdc.lattice.elasticFactor(xrdc.lattice.loadMaterials("PbTiO3"));
    a0 = xrdc.rsm.biaxialStrain(R.aPar, R.aPerp, f);
    tc.verifyEqual(R.a0, a0, 'AbsTol', 1e-6);
end

function testAnalyzeStrainRSMSyntheticPZTComposition(tc)
    % Tetragonal PZT film (a_par != a_perp) chosen so the relaxed pseudocubic
    % a0 lands on the x=0.3 Vegard a-anchor (3.99). NOTE: a *cubic* film would
    % be collinear with the cubic substrate in q-space ([1 0 3]: kPar/kPerp =
    % h/l for both lattices) and could not be separated by the auto-finder; a
    % tetragonal film gives the angular separation needed. Validates the
    % composition wiring (a0 -> interp1 on composition.a -> x).
    lambda = 1.5406; hkl = [1 0 3]; hk = hypot(hkl(1), hkl(2));
    aSub = 3.905; aParF = 3.93; aPerpF = 4.0415;   % -> a0 = 3.990 -> x = 0.3
    kSub  = [hk/aSub,  hkl(3)/aSub];
    kFilm = [hk/aParF, hkl(3)/aPerpF];
    scans = makeRsmScans(hkl, kSub, kFilm, 350, lambda);
    R = xrdc.rsm.analyzeStrainRSM(scans, Substrate="SrTiO3", Film="PZT", Reflection=hkl);
    tc.verifyEqual(R.a0, 3.99, 'AbsTol', 5e-3);
    tc.verifyEqual(R.x, 0.30, 'AbsTol', 0.02);   % interp on composition.a
    tc.verifyGreaterThan(R.relaxation, 0.15);    % partially relaxed
    tc.verifyLessThan(R.relaxation, 0.50);
end

function testAnalyzeStrainRSMRequiresAsymmetric(tc)
    scans = makeRsmScans([1 0 3], [0.256 0.768], [0.244 0.732], 300, 1.5406);
    tc.verifyError(@() xrdc.rsm.analyzeStrainRSM(scans, ...
        Substrate="SrTiO3", Film="PbTiO3", Reflection=[0 0 2]), 'xrdc:rsm:badReflection');
end

function testAnalyzeStrainRSMRealPtO2TiO2(tc)
    % GATED real-data known-answer: recover TiO2 substrate a/c from the 112
    % RSM (geometry core). Film path must run end-to-end with Film="PtO2".
    p = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data', ...
        'HP PtO2 on TiO2 001 112 RSM_C_HP PtO2 on TiO2 001 112 RSM_C.xrdml');
    tc.assumeTrue(isfile(p), 'PtO2/TiO2 112 RSM not present');
    R = xrdc.rsm.analyzeStrainRSM(p, Substrate="TiO2", Film="PtO2", Reflection=[1 1 2]);
    tc.verifyTrue(R.substrate.found);
    tc.verifyEqual(R.substrate.aMeas, 4.593, 'AbsTol', 0.02);   % rutile a
    tc.verifyEqual(R.substrate.cMeas, 2.959, 'AbsTol', 0.02);   % rutile c
    tc.verifyTrue(isfinite(R.aPar) && isfinite(R.aPerp));
end
```

- [ ] **Step 2: Run, verify failure** — `matlab -batch "addpath(pwd); addpath('tests'); run(testRsm)"`. Expected: the synthetic + error tests fail (undefined function); the real-data test is `Incomplete` if the file is absent, otherwise fails.

- [ ] **Step 3: Implement** — `+xrdc/+rsm/analyzeStrainRSM.m`:

```matlab
function R = analyzeStrainRSM(rsm, options)
%ANALYZESTRAINRSM  Strain & composition of a film from an asymmetric RSM.
%   R = xrdc.rsm.analyzeStrainRSM(rsm, Substrate=, Film=, Reflection=, ...)
%
%   Auto-detects the substrate and film peaks in an asymmetric reciprocal-
%   space map, inverts them to in-plane (a_par) and out-of-plane (a_perp)
%   lattice parameters, decomposes the film's biaxial strain to a relaxed
%   pseudocubic a0, reports the degree of relaxation against the MEASURED
%   substrate in-plane parameter, and (for an alloy film with a Vegard table,
%   e.g. PZT) the composition x. Reports numbers for the DECLARED film; it
%   does not rank candidates.
%
%   1/d convention (xrdc.rsm.toReciprocalSpace): for reflection [h k l],
%   a_par = sqrt(h^2+k^2)/kPar, a_perp = l/kPerp.
%
%   rsm        : folder/file path (-> xrdc.rsm.loadAreaScan) or a pre-loaded
%                slice struct array (each as xrdc.io.emptyScan).
%   Substrate  : declared substrate material name (default "SrTiO3").
%   Film       : declared film material name (required; must be in materials.json).
%   Reflection : 1x3 [h k l]; an asymmetric reflection (l~=0 and h^2+k^2>0).
%   NoiseFactor : forwarded to xrdc.rsm.findRsmPeaks.
%
%   Output R: see field assignments below. Caveats (pseudocubic a0; strain/
%   composition real-data-unvalidated this iteration) in SCIENTIFIC_ASSUMPTIONS.
%   Errors: xrdc:rsm:badReflection, xrdc:lattice:unknownMaterial.

    arguments
        rsm
        options.Substrate   (1,1) string  = "SrTiO3"
        options.Film        (1,1) string
        options.Reflection  (1,3) double
        options.NoiseFactor (1,1) double = 5
    end

    % Required name-value args (no default => field absent if omitted).
    if ~isfield(options, 'Reflection')
        error('xrdc:rsm:missingReflection', 'Reflection=[h k l] is required.');
    end
    if ~isfield(options, 'Film')
        error('xrdc:rsm:missingFilm', 'Film="<material>" is required.');
    end

    hkl = options.Reflection;
    hk2 = hkl(1)^2 + hkl(2)^2;
    if hkl(3) == 0 || hk2 == 0
        error('xrdc:rsm:badReflection', ...
            'Reflection [%g %g %g] is not asymmetric (need l~=0 and h^2+k^2>0).', ...
            hkl(1), hkl(2), hkl(3));
    end

    subMat  = xrdc.lattice.loadMaterials(options.Substrate);
    filmMat = xrdc.lattice.loadMaterials(options.Film);

    % --- load slices ---
    % NOTE: loadAreaScan treats a scalar string as a FOLDER and errors on a
    % file path, so a single .xrdml file must be routed through its file-list
    % branch (cell-wrapped).
    if isstruct(rsm)
        scans = rsm;
    elseif ischar(rsm) || isstring(rsm)
        rsmStr = string(rsm);
        if isfile(rsmStr)
            scans = xrdc.rsm.loadAreaScan({char(rsmStr)});
        else
            scans = xrdc.rsm.loadAreaScan(rsmStr);
        end
    else
        scans = xrdc.rsm.loadAreaScan(rsm);
    end

    % --- assemble (kPar, kPerp, intensity) cloud ---
    kp = []; kz = []; I = [];
    for j = 1:numel(scans)
        [a, b] = xrdc.rsm.toReciprocalSpace(scans(j));
        kp = [kp; a(:)]; kz = [kz; b(:)]; I = [I; double(scans(j).counts(:))]; %#ok<AGROW>
    end

    pk = xrdc.rsm.findRsmPeaks(kp, kz, I, 'NoiseFactor', options.NoiseFactor);

    sqHK = sqrt(hk2);
    invPar  = @(kpar)  sqHK / abs(kpar);
    invPerp = @(kperp) hkl(3) / kperp;

    % --- substrate ---
    subAMeas = invPar(pk.substrate.kPar);
    subCMeas = invPerp(pk.substrate.kPerp);
    % prediction from declared lattice (cubic: c=a)
    kParPred  = sqHK / subMat.a;
    kPerpPred = hkl(3) / subMat.c;
    offPct = 100 * hypot(pk.substrate.kPar - kParPred, pk.substrate.kPerp - kPerpPred) ...
                 / hypot(kParPred, kPerpPred);
    flags = pk.flags;
    if offPct > 2
        flags(end+1) = "substrateOffPrediction"; %#ok<AGROW>
    end
    R.substrate = struct('name', string(subMat.name), 'found', pk.substrate.found, ...
        'kPar', pk.substrate.kPar, 'kPerp', pk.substrate.kPerp, ...
        'aMeas', subAMeas, 'cMeas', subCMeas, 'predictOffPct', offPct);

    % --- film ---
    if ~pk.film.found
        R.film = struct('found', false, 'kPar', NaN, 'kPerp', NaN);
        [R.aPar, R.aPerp, R.cFilm, R.a0, R.strainPar, R.strainPerp, ...
         R.relaxation, R.x] = deal(NaN);
        R.pseudomorphic = false;
        R.reflection = hkl; R.elasticModel = elasticModelName(filmMat);
        R.nu = nuOf(filmMat); R.flags = flags;
        R.notes = "film peak not found; substrate-only result";
        return
    end

    aPar  = invPar(pk.film.kPar);
    aPerp = invPerp(pk.film.kPerp);
    f = xrdc.lattice.elasticFactor(filmMat);
    [a0, epsPar, epsPerp] = xrdc.rsm.biaxialStrain(aPar, aPerp, f);

    % relaxation vs MEASURED substrate in-plane parameter
    denom = a0 - subAMeas;
    if abs(denom) < 1e-6
        relax = NaN;
    else
        relax = (aPar - subAMeas) / denom;
    end
    pseudomorphic = isfinite(relax) && abs(relax) < 0.1;

    % composition (alloy films with a Vegard table only)
    x = NaN;
    if ~isempty(filmMat.composition)
        ax = filmMat.composition.a; xx = filmMat.composition.x;
        if a0 >= min(ax) && a0 <= max(ax)
            x = interp1(ax, xx, a0);
        else
            flags(end+1) = "compositionOutOfRange"; %#ok<AGROW>
        end
    end

    R.film = struct('found', true, 'kPar', pk.film.kPar, 'kPerp', pk.film.kPerp);
    R.aPar = aPar; R.aPerp = aPerp; R.cFilm = aPerp; R.a0 = a0;
    R.strainPar = epsPar; R.strainPerp = epsPerp;
    R.relaxation = relax; R.pseudomorphic = pseudomorphic; R.x = x;
    R.reflection = hkl; R.elasticModel = elasticModelName(filmMat);
    R.nu = nuOf(filmMat); R.flags = flags;
    R.notes = "a0 is a pseudocubic strain-model average (see SCIENTIFIC_ASSUMPTIONS)";
end

function name = elasticModelName(m)
    if isfield(m.elastic, 'nu') && ~isempty(m.elastic.nu)
        name = "nu";
    else
        name = "c13c33";
    end
end

function v = nuOf(m)
    if isfield(m.elastic, 'nu') && ~isempty(m.elastic.nu)
        v = m.elastic.nu;
    else
        v = NaN;     % c13/c33 model: no single nu
    end
end
```

- [ ] **Step 4: Run, verify pass** — `matlab -batch "addpath(pwd); addpath('tests'); run(testRsm)"`. Expected: synthetic + error tests pass; the real-data test passes (data is present in `data/`).

- [ ] **Step 5: Full suite + commit**

```bash
matlab -batch "r = runtests; assert(all([r.Passed]), 'suite not green')"
git add "xrdc-matlab/+xrdc/+rsm/analyzeStrainRSM.m" "xrdc-matlab/tests/testRsm.m"
git commit -m "feat(rsm): analyzeStrainRSM one-call strain/composition from an RSM"
```

---

### Task 6: Demo + documentation

**Files:**
- Create: `xrdc-matlab/examples/demoStrainRSM.m`
- Modify: `xrdc-matlab/docs/USER_GUIDE.md` (new workflow subsection)
- Modify: `xrdc-matlab/docs/SCIENTIFIC_ASSUMPTIONS.md` (new section — caveats are REQUIRED)
- Modify: `xrdc-matlab/docs/FEATURES.md` (mark done + validation-boundary note)

No new test. Verify the demo runs.

- [ ] **Step 1: Create the demo** — `examples/demoStrainRSM.m`:

```matlab
%DEMOSTRAINRSM  RSM strain & composition on the real PtO2/TiO2 112 map.
%   Run from the xrdc-matlab root. Prints the analyzeStrainRSM report.
f = fullfile('..', 'data', ...
    'HP PtO2 on TiO2 001 112 RSM_C_HP PtO2 on TiO2 001 112 RSM_C.xrdml');
if ~isfile(f)
    error('demoStrainRSM:noData', 'Place the PtO2/TiO2 112 RSM in ../data.');
end
R = xrdc.rsm.analyzeStrainRSM(f, Substrate="TiO2", Film="PtO2", Reflection=[1 1 2]);
fprintf('Substrate %s: a=%.4f c=%.4f A (found=%d, off=%.2f%%)\n', ...
    R.substrate.name, R.substrate.aMeas, R.substrate.cMeas, ...
    R.substrate.found, R.substrate.predictOffPct);
if R.film.found
    fprintf('Film: a_par=%.4f a_perp=%.4f a0=%.4f A\n', R.aPar, R.aPerp, R.a0);
    fprintf('  strain par=%+.3f%% perp=%+.3f%%  relaxation=%.2f  pseudomorphic=%d\n', ...
        100*R.strainPar, 100*R.strainPerp, R.relaxation, R.pseudomorphic);
    if ~isnan(R.x), fprintf('  composition x=%.3f\n', R.x); end
else
    fprintf('Film peak not found.\n');
end
if ~isempty(R.flags), fprintf('Flags: %s\n', strjoin(cellstr(R.flags), ', ')); end
```

- [ ] **Step 2: Run the demo** — `matlab -batch "cd('xrdc-matlab'); demoStrainRSM"`. Expected: prints TiO₂ a≈4.59 c≈2.96, a finite film line, no error.

- [ ] **Step 3: USER_GUIDE.md** — add a subsection after the existing RSM material (find the RSM/material-ID area). Insert:

```markdown
### Strain & composition from an asymmetric RSM

From a single asymmetric RSM you can recover the film's in-plane and
out-of-plane lattice parameters, its biaxial strain, the relaxed
(pseudocubic) parameter a0, the degree of relaxation, and — for an alloy
film with a Vegard table (PZT) — composition:

```matlab
R = xrdc.rsm.analyzeStrainRSM("data/your_103_RSM.xrdml", ...
        Substrate="SrTiO3", Film="PbTiO3", Reflection=[1 0 3]);
R.aPar, R.aPerp, R.a0, R.relaxation, R.x
```

`Reflection` must be asymmetric (l≠0 and h²+k²>0). The substrate peak is the
internal reference (it also absorbs a shared zero-offset). `relaxation` is 0
for a film fully strained to the substrate and 1 for a fully relaxed film.

**Trust boundary (read this):** the geometry (peak→lattice-parameter
inversion) is validated against the TiO₂ substrate of a real PtO₂/TiO₂ map.
The strain and composition outputs are, as of this release, validated only
against synthetic data — verify a reported strain or composition
independently before quoting it. For PTO/PZT, a0 is a pseudocubic average,
not a physical relaxed cubic constant; see SCIENTIFIC_ASSUMPTIONS.
```

- [ ] **Step 4: SCIENTIFIC_ASSUMPTIONS.md** — add a new section (place after the existing RSM/material-ID sections). Content must cover all six caveats from the spec:

```markdown
## RSM strain & composition (xrdc.rsm.analyzeStrainRSM)

**Convention.** Reciprocal coordinates use the 1/d convention of
`toReciprocalSpace` (|k| = 2 sinθ/λ). For reflection [h k l]:
a∥ = √(h²+k²)/kPar, a⊥ = l/kPerp.

**Biaxial decomposition.** a0 = (a⊥ + f·a∥)/(1+f),
ε∥ = (a∥−a0)/a0, ε⊥ = (a⊥−a0)/a0, with f = 2ν/(1−ν) or 2·c₁₃/c₃₃. These
satisfy ε⊥/ε∥ = −f and a = a0(1+ε) exactly.

**Validation boundary (important).** Only the geometry core (the hkl
inversion) is validated against real data — the unstrained TiO₂ substrate
peak of the PtO₂/TiO₂ 112 RSM (a=4.593, c=2.959 recovered to 4 sig figs).
Biaxial decomposition, a0 recovery, relaxation, and composition are
validated only against synthetic injected truth, which tests that the
algebra inverts the forward model, not that the model matches a real
strained film. Treat reported ε and x as unvalidated against ground truth
until a known real strained film has been run through.

**Pseudocubic approximation + Vegard consistency.** For PTO/PZT the relaxed
crystal is itself tetragonal (c/a≈1.06); a0 is a pseudocubic strain-model
average that conflates spontaneous tetragonality with epitaxial strain, not
a physical relaxed constant. Composition uses the in-plane a-axis anchors
(`composition.a`) inverted at a0 — confirm the table's a and a0 share the
same definition or x is biased. Vegard is piecewise-valid within one phase
field and breaks at the MPB.

**ν sensitivity.** a0, relaxation, and both strain components depend on ν; an
inaccurate ν biases all of them. The PtO₂ entry's ν is a placeholder.

**Tilt assumed zero.** A tilted (mosaic/miscut) film biases kPar and hence
a∥, relaxation, and composition. Not corrected this iteration.

**Centroid bias.** The analyzer/CTR streak through the substrate peak can
pull the intensity-weighted centroid; the substrate centroid is the more
exposed. Known small systematic.
```

- [ ] **Step 5: FEATURES.md** — mark RSM strain/composition as done with the validation note. Find the RSM/wanted area and add:

```markdown
- **RSM strain & composition** (`xrdc.rsm.analyzeStrainRSM`) — DONE. Auto-detects
  substrate+film peaks in an asymmetric RSM, recovers a∥/a⊥, biaxial strain,
  relaxed a0, relaxation, and PZT composition for a declared film. Geometry
  validated on real TiO₂ substrate; strain/composition synthetic-only this
  iteration (see SCIENTIFIC_ASSUMPTIONS). Candidate ranking / identifyMaterial
  integration is a deliberate follow-on.
```

- [ ] **Step 6: Full suite + commit**

```bash
matlab -batch "r = runtests; assert(all([r.Passed]), 'suite not green')"
git add "xrdc-matlab/examples/demoStrainRSM.m" "xrdc-matlab/docs/USER_GUIDE.md" "xrdc-matlab/docs/SCIENTIFIC_ASSUMPTIONS.md" "xrdc-matlab/docs/FEATURES.md"
git commit -m "docs(rsm): demo + USER_GUIDE/SCIENTIFIC_ASSUMPTIONS/FEATURES for strain RSM"
```

---

## Notes for the executor

- **Order matters:** Task 1 (elasticFactor) and Task 3 (PtO₂) are prerequisites for Task 5. Task 2 and Task 4 are independent of each other but both feed Task 5. Keep the plan order.
- **Real-data test is NOT gated-away here:** the PtO₂/TiO₂ 112 RSM is present in `data/`, so `testAnalyzeStrainRSMRealPtO2TiO2` must actually pass (substrate a/c within 0.02 Å). If it only reports `Incomplete`, the file path is wrong — fix it, don't accept the skip.
- **The strain label trap (Task 2):** the implementation is the compact form, which is correct by construction. Do not "expand" it into the explicit ±(1−ν)/(1+ν) form during cleanup; if you do, the sign test in Task 2 is what guards you.
- **Do not** wire `analyzeStrainRSM` into `identifyMaterial` or add candidate ranking — explicitly out of scope this iteration.
- Report DONE / DONE_WITH_CONCERNS / BLOCKED per task.
```
