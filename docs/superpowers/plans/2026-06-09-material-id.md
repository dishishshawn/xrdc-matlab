# Material Identification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Identify which material(s) (PZT, SrTiO3, PbTiO3, SrRuO3) produced the (00l) peaks of a θ-2θ scan, with strain-aware matching, and report c, strain, relaxation, confidence, and (for PZT) a strain-corrected composition estimate.

**Architecture:** Peaks → Kβ/W-Lα ghost filter → claim declared-substrate series → cluster remaining peaks into harmonic (00l) series → least-squares c per series → score each series' c against each candidate material's [c_bulk, c_pred(pseudomorphic)] interval → ranked candidate sets. Spec: `docs/superpowers/specs/2026-06-09-material-id-design.md` — read it first.

**Tech Stack:** MATLAB R2022b+ (`arguments` blocks, `string`), `functiontests` unit tests, no new toolbox dependencies. All work happens under the repo's `xrdc-matlab/` subfolder (note: paths below are relative to repo root `C:\Users\TNTMi\Desktop\xrdc-matlab`).

**Conventions that apply to every task** (from `xrdc-matlab/CLAUDE.md`):
- Error IDs: `xrdc:<subpkg>:<reason>`; warnings likewise.
- Angles in degrees in public APIs, wavelengths/lattice in Å.
- One function per file; `camelCase`; docstring header in the house style (see `+xrdc/+lattice/twoThetaToD.m` for the template).
- Run the full suite (`runtests` from `xrdc-matlab/`) before each commit; commit but **never push**.

**File map:**

| File | Role |
|---|---|
| `xrdc-matlab/+xrdc/+data/materials.json` | Create — material database |
| `xrdc-matlab/+xrdc/+lattice/loadMaterials.m` | Create — load/normalize/lookup database |
| `xrdc-matlab/+xrdc/+peaks/filterGhostPeaks.m` | Create — Cu Kβ / W Lα ghost filter |
| `xrdc-matlab/+xrdc/+lattice/groupHarmonicSeries.m` | Create — cluster d-spacings into (00l) series |
| `xrdc-matlab/+xrdc/+lattice/identifyMaterial.m` | Create — full pipeline |
| `xrdc-matlab/tests/testIdentify.m` | Create — all tests for this feature (feature-focused suite; spans +peaks and +lattice deliberately) |
| `xrdc-matlab/build/buildStandalone.m` | Modify line 96 — bundle materials.json |
| `xrdc-matlab/docs/FEATURES.md` | Modify — feature list entry |
| `xrdc-matlab/xrdcApp.m` | Modify — GUI wiring |

---

### Task 1: Material database + loader

**Files:**
- Create: `xrdc-matlab/+xrdc/+data/materials.json`
- Create: `xrdc-matlab/+xrdc/+lattice/loadMaterials.m`
- Create: `xrdc-matlab/tests/testIdentify.m`

- [ ] **Step 1: Write the database**

Every entry has the SAME top-level fields (jsondecode only forms a struct array when field names match across all array elements — that's why `composition: null` appears on non-PZT entries). Lattice values: STO from `substrates.json` (Substrates.def lineage); PTO/PZT tetragonal-side values are the committed v1 anchors per the spec (approximate ±0.01 Å — adequate for ID windows, not precision composition; see refs strings). ν values are the standard biaxial-strain Poisson ratios used in perovskite film analyses; the code also supports `{"c13":..,"c33":..}` for entries that later get measured cij.

```json
{
  "_comment": "Material reference data for xrdc.lattice.identifyMaterial. Lattice values in Angstrom. role: substrate|film|both. elastic: either {nu} (cubic/biaxial approx, factor 2*nu/(1-nu)) or {c13,c33} in GPa (tetragonal, factor 2*c13/c33). composition (PZT only): tetragonal-side anchors for a(x), c(x), x = Zr fraction. All entries share identical top-level fields so jsondecode yields a struct array.",
  "materials": [
    {
      "name": "SrTiO3",
      "aliases": ["STO", "strontium titanate"],
      "system": "cubic",
      "a": 3.905,
      "c": 3.905,
      "role": "both",
      "elastic": { "nu": 0.232 },
      "composition": null,
      "refs": "a = 3.905 A (Substrates.def / substrates.json); nu = 0.232 (commonly used STO biaxial Poisson ratio)"
    },
    {
      "name": "SrRuO3",
      "aliases": ["SRO", "strontium ruthenate"],
      "system": "cubic",
      "a": 3.93,
      "c": 3.93,
      "role": "film",
      "elastic": { "nu": 0.30 },
      "composition": null,
      "refs": "pseudocubic a = 3.93 A (orthorhombic SRO treated pseudocubically, standard for (00l) film work); nu = 0.30 nominal perovskite value"
    },
    {
      "name": "PbTiO3",
      "aliases": ["PTO", "lead titanate"],
      "system": "tetragonal",
      "a": 3.904,
      "c": 4.152,
      "role": "film",
      "elastic": { "nu": 0.30 },
      "composition": null,
      "refs": "a = 3.904, c = 4.152 A (Shirane & Hoshino 1951 tetragonal PbTiO3); nu = 0.30 nominal"
    },
    {
      "name": "PZT",
      "aliases": ["Pb(Zr,Ti)O3", "PbZrTiO3", "lead zirconate titanate"],
      "system": "tetragonal",
      "a": 4.036,
      "c": 4.146,
      "role": "film",
      "elastic": { "nu": 0.30 },
      "composition": {
        "x": [0.0, 0.1, 0.2, 0.3, 0.4, 0.52],
        "a": [3.904, 3.93, 3.96, 3.99, 4.01, 4.036],
        "c": [4.152, 4.15, 4.149, 4.148, 4.147, 4.146]
      },
      "refs": "nominal entry = MPB 52/48; composition anchors interpolated from the tetragonal-side phase data of Jaffe, Cook & Jaffe, Piezoelectric Ceramics (1971) / Shirane & Suzuki (1952); approximate +-0.01 A, adequate for ID windows"
    }
  ]
}
```

- [ ] **Step 2: Write the failing tests**

Create `xrdc-matlab/tests/testIdentify.m`:

```matlab
function tests = testIdentify
%TESTIDENTIFY  Tests for material identification (materials.json,
%   filterGhostPeaks, groupHarmonicSeries, identifyMaterial).
    tests = functiontests(localfunctions);
end

% ---------- loadMaterials ----------

function testLoadMaterialsAll(tc)
    M = xrdc.lattice.loadMaterials();
    tc.verifyEqual(numel(M), 4);
    tc.verifyTrue(all(isfield(M, {'name','aliases','system','a','c','role','elastic','composition'})));
    names = string({M.name});
    tc.verifyTrue(all(ismember(["SrTiO3","SrRuO3","PbTiO3","PZT"], names)));
end

function testLoadMaterialsLookupByAlias(tc)
    e = xrdc.lattice.loadMaterials("PTO");
    tc.verifyEqual(string(e.name), "PbTiO3");
    tc.verifyEqual(e.c, 4.152, 'AbsTol', 1e-9);
    e2 = xrdc.lattice.loadMaterials("srtio3");   % case-insensitive
    tc.verifyEqual(e2.a, 3.905, 'AbsTol', 1e-9);
end

function testLoadMaterialsUnknownErrors(tc)
    tc.verifyError(@() xrdc.lattice.loadMaterials("kryptonite"), ...
        'xrdc:lattice:unknownMaterial');
end
```

- [ ] **Step 3: Run tests to verify they fail**

From `xrdc-matlab/`: `matlab -batch "results = runtests('tests/testIdentify.m'); disp(table(results)); assertSuccess(results)"`
Expected: FAIL — `xrdc.lattice.loadMaterials` undefined.

- [ ] **Step 4: Implement the loader**

`xrdc-matlab/+xrdc/+lattice/loadMaterials.m`:

```matlab
function M = loadMaterials(name)
%LOADMATERIALS  Load the material database (+xrdc/+data/materials.json).
%   M = xrdc.lattice.loadMaterials()        — struct array, all entries
%   e = xrdc.lattice.loadMaterials(name)    — single entry by name/alias
%
%   Lookup is case-insensitive over .name and .aliases. Errors with
%   xrdc:lattice:unknownMaterial (message lists valid names) on no match.
%
%   See also xrdc.lattice.identifyMaterial.

    arguments
        name (1,1) string = ""
    end

    persistent cache
    if isempty(cache)
        jsonPath = fullfile(fileparts(mfilename('fullpath')), ...
            '..', '+data', 'materials.json');
        raw = jsondecode(fileread(jsonPath));
        mats = raw.materials;
        if iscell(mats)   % defensive: heterogeneous fields decode as cell
            mats = [mats{:}];
        end
        cache = mats(:).';
    end

    if name == ""
        M = cache;
        return
    end

    target = lower(name);
    for e = cache
        candidates = lower([string(e.name); string(e.aliases(:))]);
        if any(candidates == target)
            M = e;
            return
        end
    end
    error('xrdc:lattice:unknownMaterial', ...
        'Unknown material "%s". Valid names: %s.', name, ...
        strjoin(string({cache.name}), ', '));
end
```

- [ ] **Step 5: Run tests to verify they pass**

Same command as Step 3. Expected: the three loadMaterials tests PASS.

- [ ] **Step 6: Commit**

```bash
git add xrdc-matlab/+xrdc/+data/materials.json xrdc-matlab/+xrdc/+lattice/loadMaterials.m xrdc-matlab/tests/testIdentify.m
git commit -m "feat: material database + loader for material identification"
```

---

### Task 2: Kβ / W-Lα ghost filter

**Files:**
- Create: `xrdc-matlab/+xrdc/+peaks/filterGhostPeaks.m`
- Modify: `xrdc-matlab/tests/testIdentify.m` (append)

Physics: a strong reflection with spacing d (computed from the main λ, Cu Kα1 1.5406 Å) also diffracts the tube's Cu Kβ (1.3922 Å) and tungsten-contamination W Lα1 (1.4763 Å) lines, producing small "ghost" peaks at 2θ_ghost = 2·asin(λ_ghost/2d) — always BELOW the parent peak. The wavelengths correspond to xrayLines.json energies (Cu-Kbeta1 8905.413 eV, W-Lalpha1 8398.242 eV).

- [ ] **Step 1: Write the failing tests** (append to `tests/testIdentify.m`)

```matlab
% ---------- filterGhostPeaks ----------

function testGhostFilterRemovesKbeta(tc)
    lambda = 1.5406;
    % STO (002): d = 3.905/2; Kbeta ghost of that reflection:
    d002   = 3.905/2;
    parent = xrdc.lattice.dToTwoTheta(d002, lambda);     % ~46.48
    ghost  = xrdc.lattice.dToTwoTheta(d002, 1.3922);     % ~41.77
    tt = [ghost; parent];  counts = [5e4; 1e6];          % ghost is 5% of parent
    [keep, ghosts] = xrdc.peaks.filterGhostPeaks(tt, counts, lambda);
    tc.verifyEqual(keep, [false; true]);
    tc.verifyEqual(height(ghosts), 1);
    tc.verifyEqual(ghosts.twoTheta(1), ghost, 'AbsTol', 1e-9);
end

function testGhostFilterRespectsIntensityRatio(tc)
    % A peak AT the ghost position but as strong as the parent is real.
    lambda = 1.5406;
    d002   = 3.905/2;
    tt = [xrdc.lattice.dToTwoTheta(d002, 1.3922); ...
          xrdc.lattice.dToTwoTheta(d002, lambda)];
    counts = [9e5; 1e6];   % 90% of parent: not a ghost
    keep = xrdc.peaks.filterGhostPeaks(tt, counts, lambda);
    tc.verifyEqual(keep, [true; true]);
end

function testGhostFilterNoIntensitiesWarns(tc)
    tt = [41.77; 46.48];
    f = @() xrdc.peaks.filterGhostPeaks(tt, nan(2,1), 1.5406);
    tc.verifyWarning(f, 'xrdc:peaks:noIntensity');
    keep = f();  %#ok<NASGU>  % suppress duplicate warning display
    [keep, ~] = xrdc.peaks.filterGhostPeaks(tt, nan(2,1), 1.5406);
    tc.verifyEqual(keep, [true; true]);   % nothing removed blind
end
```

- [ ] **Step 2: Run tests to verify they fail**

`matlab -batch "results = runtests('tests/testIdentify.m'); disp(table(results)); assertSuccess(results)"`
Expected: FAIL — `xrdc.peaks.filterGhostPeaks` undefined.

- [ ] **Step 3: Implement**

`xrdc-matlab/+xrdc/+peaks/filterGhostPeaks.m`:

```matlab
function [keep, ghosts] = filterGhostPeaks(twoTheta, counts, lambda, opts)
%FILTERGHOSTPEAKS  Flag Cu-Kbeta / W-Lalpha ghost peaks of strong reflections.
%   [keep, ghosts] = xrdc.peaks.filterGhostPeaks(twoTheta, counts, lambda)
%
%   A reflection with spacing d (from the main wavelength) also diffracts
%   the Cu Kbeta (1.3922 A) and W Lalpha1 (1.4763 A) tube lines, leaving
%   small satellite peaks at 2*asin(lambdaGhost/2d), below the parent.
%   Any peak within PositionTol of a predicted ghost position whose
%   intensity is at most MaxRatio x the parent is flagged.
%
%   Inputs
%     twoTheta : peak positions, degrees (vector)
%     counts   : peak intensities (same size; all-NaN -> warning
%                xrdc:peaks:noIntensity and nothing is flagged, because
%                strong/weak cannot be ranked)
%     lambda   : main wavelength in A
%
%   Name-Value
%     GhostLambdas (default [1.3922 1.4763])  — Cu Kbeta1, W Lalpha1 (A)
%     PositionTol  (default 0.15)             — degrees 2-theta
%     MaxRatio     (default 0.3)              — ghost/parent intensity cap
%
%   Outputs
%     keep   : logical mask, true = real peak
%     ghosts : table (twoTheta, counts, parentTwoTheta, ghostLambda)
%
%   See also xrdc.lattice.identifyMaterial.

    arguments
        twoTheta (:,1) double
        counts   (:,1) double
        lambda   (1,1) double {mustBePositive}
        opts.GhostLambdas (1,:) double = [1.3922, 1.4763]
        opts.PositionTol  (1,1) double {mustBePositive} = 0.15
        opts.MaxRatio     (1,1) double {mustBePositive} = 0.3
    end

    n = numel(twoTheta);
    assert(numel(counts) == n, 'xrdc:peaks:sizeMismatch', ...
        'twoTheta and counts must be the same length.');
    keep = true(n, 1);
    ghosts = table('Size', [0 4], ...
        'VariableTypes', {'double','double','double','double'}, ...
        'VariableNames', {'twoTheta','counts','parentTwoTheta','ghostLambda'});

    if all(isnan(counts))
        warning('xrdc:peaks:noIntensity', ...
            ['No peak intensities provided - ghost filtering skipped ', ...
             '(cannot rank strong vs weak peaks).']);
        return
    end

    % Strongest first: a ghost can never claim its own parent.
    [~, order] = sort(counts, 'descend', 'MissingPlacement', 'last');
    for p = order(:).'
        if ~keep(p) || isnan(counts(p)), continue, end
        d = xrdc.lattice.twoThetaToD(twoTheta(p), lambda);
        for gl = opts.GhostLambdas
            g2t = xrdc.lattice.dToTwoTheta(d, gl);
            if isnan(g2t), continue, end
            hits = find(keep & (1:n).' ~= p ...
                & abs(twoTheta - g2t) <= opts.PositionTol ...
                & counts <= opts.MaxRatio * counts(p));
            for j = hits(:).'
                keep(j) = false;
                ghosts(end+1, :) = {twoTheta(j), counts(j), twoTheta(p), gl}; %#ok<AGROW>
            end
        end
    end
end
```

- [ ] **Step 4: Run tests to verify they pass**

Same command. Expected: all testIdentify tests so far PASS.

- [ ] **Step 5: Commit**

```bash
git add xrdc-matlab/+xrdc/+peaks/filterGhostPeaks.m xrdc-matlab/tests/testIdentify.m
git commit -m "feat: Cu-Kbeta / W-Lalpha ghost peak filter"
```

---

### Task 3: Harmonic series grouping

**Files:**
- Create: `xrdc-matlab/+xrdc/+lattice/groupHarmonicSeries.m`
- Modify: `xrdc-matlab/tests/testIdentify.m` (append)

Each (00l) family produces d(00l) = c/l. Greedy clustering: seed a candidate c from every (peak, assumed-order) pair, collect peaks consistent with that c, take the biggest/tightest series, repeat. Note the greedy pass naturally prefers the SMALLEST c that explains a set (a {002,004} pair groups as orders {1,2} of c/2) — the order-doubling correction happens at the database-matching stage in Task 4, where plausibility can decide.

- [ ] **Step 1: Write the failing tests** (append to `tests/testIdentify.m`)

```matlab
% ---------- groupHarmonicSeries ----------

function testGroupTwoCleanSeries(tc)
    cA = 3.905;  cB = 4.151;
    d  = [cA./(1:4), cB./(1:3)].';
    [S, singles] = xrdc.lattice.groupHarmonicSeries(d);
    tc.verifyEqual(numel(S), 2);
    tc.verifyEmpty(singles);
    cs = sort([S.c]);
    tc.verifyEqual(cs, sort([cA cB]), 'AbsTol', 1e-9);
    tc.verifyEqual(sort(S([S.c] == cA).orders), 1:4);
end

function testGroupLeavesSingles(tc)
    d = [3.905./(1:3), 2.31].';   % 2.31 A fits no order of 3.905 within tol
    [S, singles] = xrdc.lattice.groupHarmonicSeries(d);
    tc.verifyEqual(numel(S), 1);
    tc.verifyEqual(singles, 4);
end

function testGroupToleranceRejects(tc)
    % Second "order" 1.2% off — outside the 0.5% default tolerance.
    d = [3.905; 3.905/2 * 1.012];
    [S, singles] = xrdc.lattice.groupHarmonicSeries(d);
    tc.verifyEmpty(S);
    tc.verifyEqual(sort(singles(:).'), [1 2]);
end
```

- [ ] **Step 2: Run tests to verify they fail**

Same runtests command. Expected: FAIL — `groupHarmonicSeries` undefined.

- [ ] **Step 3: Implement**

`xrdc-matlab/+xrdc/+lattice/groupHarmonicSeries.m`:

```matlab
function [S, singles] = groupHarmonicSeries(d, opts)
%GROUPHARMONICSERIES  Cluster d-spacings into harmonic (00l) series.
%   [S, singles] = xrdc.lattice.groupHarmonicSeries(d)
%
%   In a (00l)-oriented scan every phase gives d(00l) = c/l. Greedy
%   clustering: every (peak, assumed order 1..MaxOrder) pair seeds a
%   candidate c; peaks whose nearest-integer order reproduces their d
%   within Tolerance join; the largest (then tightest) series is taken
%   and the loop repeats on the remainder. One peak per order (closest
%   wins). The greedy pass prefers the smallest c that explains a set —
%   order-doubling disambiguation belongs to the database-matching stage
%   (see identifyMaterial).
%
%   Inputs
%     d : vector of d-spacings (A)
%
%   Name-Value
%     Tolerance (default 0.005) — max relative |d - c/l| / (c/l)
%     MaxOrder  (default 4)
%     CRange    (default [2 8]) — plausible c window (A) for seeds
%
%   Outputs
%     S       : struct array, fields .members (indices into d), .orders,
%               .c (least-squares), .cSigma, .residRms (relative),
%               .evenOnly (all matched orders even)
%     singles : indices of peaks in no multi-peak series
%
%   See also xrdc.lattice.identifyMaterial.

    arguments
        d (:,1) double {mustBePositive}
        opts.Tolerance (1,1) double {mustBePositive} = 0.005
        opts.MaxOrder  (1,1) double {mustBeInteger, mustBePositive} = 4
        opts.CRange    (1,2) double = [2 8]
    end

    n = numel(d);
    unassigned = true(n, 1);
    S = struct('members', {}, 'orders', {}, 'c', {}, 'cSigma', {}, ...
               'residRms', {}, 'evenOnly', {});

    while true
        best = struct('members', [], 'orders', [], 'rms', Inf);
        idx = find(unassigned).';
        for j = idx
            for l = 1:opts.MaxOrder
                c0 = d(j) * l;
                if c0 < opts.CRange(1) || c0 > opts.CRange(2), continue, end
                [mem, ord, rms] = collectMembers(d, idx, c0, opts);
                better = numel(mem) > numel(best.members) || ...
                        (numel(mem) == numel(best.members) && rms < best.rms);
                if numel(mem) >= 2 && better
                    best = struct('members', mem, 'orders', ord, 'rms', rms);
                end
            end
        end
        if numel(best.members) < 2, break, end

        cVals  = d(best.members) .* best.orders(:);
        c      = mean(cVals);
        sigma  = std(cVals) / sqrt(numel(cVals));
        resid  = (d(best.members) - c ./ best.orders(:)) ./ (c ./ best.orders(:));
        S(end+1) = struct('members', best.members, 'orders', best.orders, ...
            'c', c, 'cSigma', sigma, 'residRms', sqrt(mean(resid.^2)), ...
            'evenOnly', all(mod(best.orders, 2) == 0)); %#ok<AGROW>
        unassigned(best.members) = false;
    end

    singles = find(unassigned);
end

function [mem, ord, rms] = collectMembers(d, idx, c0, opts)
%COLLECTMEMBERS  Peaks consistent with series constant c0; one per order.
    mem = []; ord = []; resid = [];
    for k = idx
        lk = round(c0 / d(k));
        if lk < 1 || lk > opts.MaxOrder, continue, end
        dPred = c0 / lk;
        r = abs(d(k) - dPred) / dPred;
        if r > opts.Tolerance, continue, end
        prev = find(ord == lk, 1);
        if isempty(prev)
            mem(end+1) = k; ord(end+1) = lk; resid(end+1) = r; %#ok<AGROW>
        elseif r < resid(prev)          % closer claimant wins the order
            mem(prev) = k; resid(prev) = r;
        end
    end
    rms = sqrt(mean(resid.^2));
    if isempty(mem), rms = Inf; end
end
```

- [ ] **Step 4: Run tests to verify they pass**

Same command. Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add xrdc-matlab/+xrdc/+lattice/groupHarmonicSeries.m xrdc-matlab/tests/testIdentify.m
git commit -m "feat: harmonic (00l) series clustering"
```

---

### Task 4: identifyMaterial pipeline

**Files:**
- Create: `xrdc-matlab/+xrdc/+lattice/identifyMaterial.m`
- Modify: `xrdc-matlab/tests/testIdentify.m` (append)

The strain model, per spec: ε∥ = (a_sub − a_bulk)/a_bulk; ε⊥ = −f·ε∥ with f = 2ν/(1−ν) (elastic has `nu`) or f = 2·c13/c33 (elastic has `c13`,`c33`); c_pred = c_bulk·(1+ε⊥). Score c_meas against the hull [min,max] of {c_bulk, c_pred} (for PZT: hull over the whole composition grid), padded by WindowPad. **Never a lone winner:** all accepted candidates are returned ranked. Order ambiguity: every series is also evaluated as 2c (orders doubled, when 2·max(orders) ≤ MaxOrder) and c/2 (when all orders even); the hypothesis with the best-scoring candidate wins, flags record what happened.

Worked numbers used by the tests (verify these before writing assertions — λ = 1.5406 Å):
- STO (001..004): 2θ = 22.756, 46.476, 72.573, 104.175°
- PTO pseudomorphic on STO: ε∥ = +2.56e-4, f = 2·0.3/0.7 = 0.857, c_pred = 4.1511 Å → (001..004) at 21.385, 43.567, 67.654, 95.851°
- PZT 52/48 pseudomorphic on STO: ε∥ = −3.246%, c_pred = 4.146·1.0278 = 4.2613 Å
- Ranking tie-breaks: score desc, then |c_meas − c_pred| asc, then fixed-composition entries before composition-model entries (parsimony — dilute PZT can imitate PTO arbitrarily well, so PTO must win the tie and the ambiguity flag must be set).

- [ ] **Step 1: Write the failing tests** (append to `tests/testIdentify.m`)

```matlab
% ---------- identifyMaterial: synthetic known-answer ----------

function tt = stoPlusFilmPeaks(cFilm, ordersFilm)
%Helper (not a test): STO substrate (001..004) + film series 2-thetas.
    lambda = 1.5406;
    ttSub  = xrdc.lattice.dToTwoTheta(3.905 ./ (1:4), lambda);
    ttFilm = xrdc.lattice.dToTwoTheta(cFilm ./ ordersFilm, lambda);
    tt = [ttSub(:); ttFilm(:)];
end

function testIdentifyPseudomorphicPTO(tc)
    tt = stoPlusFilmPeaks(4.1511, 1:3);
    R  = xrdc.lattice.identifyMaterial(tt, 1.5406, Substrate="SrTiO3");
    tc.verifyTrue(R.substrate.found);
    tc.verifyEqual(numel(R.substrate.twoTheta), 4);
    tc.verifyEqual(height(R.series), 1);
    tc.verifyEqual(R.series.bestMatch(1), "PbTiO3");
    tc.verifyEqual(R.series.cMeas(1), 4.1511, 'AbsTol', 2e-3);
    % Dilute PZT imitates PTO: both accepted, ambiguity flagged, PTO first.
    cand = R.series.candidates{1};
    tc.verifyTrue(any(string({cand.name}) == "PZT"));
    tc.verifyTrue(any(R.series.flags{1} == "ambiguous"));
    % Pseudomorphic: relaxation near 0 for the PTO candidate.
    pto = cand(string({cand.name}) == "PbTiO3");
    tc.verifyEqual(pto.relaxation, 0, 'AbsTol', 0.15);
end

function testIdentifyPseudomorphicPZTUnambiguous(tc)
    % c = 4.2613 only PZT can explain (PTO window ends ~4.152 + pad).
    tt = stoPlusFilmPeaks(4.2613, 1:3);
    R  = xrdc.lattice.identifyMaterial(tt, 1.5406);
    tc.verifyEqual(R.series.bestMatch(1), "PZT");
    cand = R.series.candidates{1};
    tc.verifyFalse(any(string({cand.name}) == "PbTiO3"));
    pzt = cand(string({cand.name}) == "PZT");
    tc.verifyEqual(pzt.x, 0.52, 'AbsTol', 0.06);   % strain-corrected estimate
    tc.verifyFalse(any(R.series.flags{1} == "ambiguous"));
end

function testIdentifyRelaxedPZTNoComposition(tc)
    % Fully relaxed MPB PZT: c = bulk 4.146; pseudomorphic inversion must
    % refuse (x = NaN) because c_meas is below the pseudomorphic range.
    tt = stoPlusFilmPeaks(4.146, 1:3);
    R  = xrdc.lattice.identifyMaterial(tt, 1.5406);
    cand = R.series.candidates{1};
    pzt  = cand(string({cand.name}) == "PZT");
    tc.verifyTrue(isnan(pzt.x));
    tc.verifyEqual(pzt.relaxation, 1, 'AbsTol', 0.2);
end

function testIdentifyOrderDoubling(tc)
    % Film shows only 002 and 004 (001/003 too weak): greedy grouping sees
    % orders {1,2} of c/2 = 2.0756 — database matching must recover c.
    tt = stoPlusFilmPeaks(4.1511, [2 4]);
    R  = xrdc.lattice.identifyMaterial(tt, 1.5406);
    tc.verifyEqual(R.series.bestMatch(1), "PbTiO3");
    tc.verifyEqual(R.series.cMeas(1), 4.1511, 'AbsTol', 2e-3);
    tc.verifyTrue(any(R.series.flags{1} == "orderDoubled"));
end

function testIdentifyGhostNotASeries(tc)
    % Kbeta ghost of STO 002 injected with intensities: must be filtered,
    % not grouped or named.
    lambda = 1.5406;
    tt = stoPlusFilmPeaks(4.1511, 1:3);
    ghost2t = xrdc.lattice.dToTwoTheta(3.905/2, 1.3922);
    pk = struct('twoTheta', num2cell([tt; ghost2t]), ...
                'counts',   num2cell([1e5*ones(numel(tt),1); 2e4]));
    % make the parent (STO 002, element 2) strong:
    pk(2).counts = 1e6;
    R = xrdc.lattice.identifyMaterial(pk, lambda);
    tc.verifyEqual(height(R.ghosts), 1);
    tc.verifyEqual(height(R.series), 1);            % just the PTO series
    tc.verifyEqual(R.series.bestMatch(1), "PbTiO3");
end

function testIdentifySubstrateMissingWarns(tc)
    ttFilmOnly = xrdc.lattice.dToTwoTheta(4.1511 ./ (1:3), 1.5406);
    tc.verifyWarning(@() xrdc.lattice.identifyMaterial(ttFilmOnly(:), 1.5406), ...
        'xrdc:lattice:substrateNotFound');
end

function testIdentifyErrorPaths(tc)
    tc.verifyError(@() xrdc.lattice.identifyMaterial([], 1.5406), ...
        'xrdc:lattice:noPeaks');
    tc.verifyError(@() xrdc.lattice.identifyMaterial([22.7; 46.5], 1.5406, ...
        Substrate="kryptonite"), 'xrdc:lattice:unknownMaterial');
end
```

NOTE: the helper `stoPlusFilmPeaks` is a local function but NOT a test (functiontests only treats `test*` functions as tests; lowercase helper is fine).

- [ ] **Step 2: Run tests to verify they fail**

Same command. Expected: FAIL — `identifyMaterial` undefined.

- [ ] **Step 3: Implement**

`xrdc-matlab/+xrdc/+lattice/identifyMaterial.m`:

```matlab
function R = identifyMaterial(peaks, lambda, opts)
%IDENTIFYMATERIAL  Identify material(s) from (00l) peak positions.
%   R = xrdc.lattice.identifyMaterial(peaks, lambda, Name=Value)
%
%   Pipeline (see docs/superpowers/specs/2026-06-09-material-id-design.md):
%   ghost filter -> claim declared substrate series -> harmonic grouping
%   -> per-series least-squares c -> strain-aware candidate matching.
%   Candidates are scored against the interval [c_bulk, c_pred] where
%   c_pred is the pseudomorphic-on-substrate prediction
%   (eps_perp = -f*eps_par, f = 2nu/(1-nu) or 2*c13/c33); a film sits
%   anywhere on that interval depending on relaxation. ALL accepted
%   candidates are returned ranked - never a silent single winner.
%
%   Inputs
%     peaks  : one of
%              - vector of peak 2-theta positions (degrees)
%              - struct array from xrdc.peaks.findPeaks (.twoTheta, .counts)
%              - scan struct (vector .twoTheta/.counts) -> findPeaks is run
%     lambda : wavelength in A (default Cu Kalpha1 1.5406)
%
%   Name-Value
%     Substrate    (default "SrTiO3") - declared substrate (name or alias)
%     SubstrateTol (default 0.001)    - relative d window to claim substrate
%     Tolerance    (default 0.005)    - harmonic-series grouping tolerance
%     MaxOrder     (default 4)
%     WindowPad    (default 0.015)    - candidate window pad (relative c)
%     GhostFilter  (default true)
%
%   Output R (struct):
%     .substrate : .name, .found, .twoTheta, .cMeas, .cRef
%     .series    : table - cMeas, cSigma, orders, twoTheta, bestMatch,
%                  bestScore, candidates (cell of struct arrays with
%                  .name .score .misfit .cBulk .cPred .strainVsBulk
%                  .relaxation .x), flags (cell of string arrays)
%     .unassigned: 2-theta of peaks in no series and matching nothing
%     .ghosts    : table from xrdc.peaks.filterGhostPeaks
%     .notes     : string array (incl. the PZT strain-confounding caveat)
%
%   See also xrdc.lattice.loadMaterials, xrdc.lattice.groupHarmonicSeries,
%            xrdc.peaks.filterGhostPeaks, xrdc.peaks.findPeaks.

    arguments
        peaks
        lambda (1,1) double {mustBePositive} = 1.5406
        opts.Substrate    (1,1) string = "SrTiO3"
        opts.SubstrateTol (1,1) double {mustBePositive} = 0.001
        opts.Tolerance    (1,1) double {mustBePositive} = 0.005
        opts.MaxOrder     (1,1) double {mustBeInteger, mustBePositive} = 4
        opts.WindowPad    (1,1) double {mustBePositive} = 0.015
        opts.GhostFilter  (1,1) logical = true
    end

    PZT_CAVEAT = "Composition from (00l) alone is strain-confounded; the x" + ...
        " estimate assumes a fully pseudomorphic film. Deconvolving strain" + ...
        " and composition needs the in-plane parameter (RSM or asymmetric" + ...
        " reflection).";

    [tt, counts] = normalizeInput(peaks);
    if isempty(tt)
        error('xrdc:lattice:noPeaks', 'No peaks supplied or detected.');
    end

    % --- 1. ghost filter -------------------------------------------------
    ghosts = table();
    if opts.GhostFilter
        [keep, ghosts] = xrdc.peaks.filterGhostPeaks(tt, counts, lambda);
        tt = tt(keep); counts = counts(keep); %#ok<NASGU>
    end

    % --- 2. substrate: declared input, confirmed not discovered ----------
    sub = xrdc.lattice.loadMaterials(opts.Substrate);
    if ~ismember(string(sub.role), ["substrate", "both"])
        error('xrdc:lattice:badSubstrate', ...
            '%s has role "%s" - not usable as a substrate.', sub.name, sub.role);
    end
    d = xrdc.lattice.twoThetaToD(tt, lambda);
    cSub = sub.c;
    isSubPeak = false(size(d));
    subOrders = [];
    for l = 1:opts.MaxOrder
        dPred = cSub / l;
        [err, j] = min(abs(d - dPred) / dPred);
        if err <= opts.SubstrateTol && ~isSubPeak(j)
            isSubPeak(j) = true;
            subOrders(end+1) = l; %#ok<AGROW>
        end
    end
    R.substrate = struct('name', string(sub.name), 'found', any(isSubPeak), ...
        'twoTheta', tt(isSubPeak), 'cRef', cSub, ...
        'cMeas', meanOrNaN(d(isSubPeak) .* subOrders(:)));
    if ~R.substrate.found
        warning('xrdc:lattice:substrateNotFound', ...
            ['Declared substrate %s: no (00l) series found within %.2f%%. ', ...
             'Proceeding with all peaks unassigned.'], ...
            sub.name, 100 * opts.SubstrateTol);
    end

    % --- 3. group remaining peaks into harmonic series --------------------
    ttF = tt(~isSubPeak);  dF = d(~isSubPeak);
    [S, singleIdx] = xrdc.lattice.groupHarmonicSeries(dF, ...
        Tolerance=opts.Tolerance, MaxOrder=opts.MaxOrder);

    % --- 4. strain-aware naming, with order-ambiguity hypotheses ----------
    M = xrdc.lattice.loadMaterials();
    films = M(ismember(string({M.role}), ["film", "both"]));
    rows = cell(0, 8);
    unassigned = [];
    for s = S
        hyps = orderHypotheses(s, opts.MaxOrder);
        bestH = []; bestCand = []; bestScoreH = -Inf; flags = strings(0,1);
        anyAcceptedCount = 0;
        for h = hyps
            cand = rankCandidates(films, h.c, sub.a, opts.WindowPad);
            if ~isempty(cand), anyAcceptedCount = anyAcceptedCount + 1; end
            hScore = -Inf;
            if ~isempty(cand), hScore = cand(1).score; end
            if hScore > bestScoreH || isempty(bestH)
                bestH = h; bestCand = cand; bestScoreH = hScore;
            end
        end
        if bestH.kind ~= "asIs", flags(end+1) = bestH.kind; end %#ok<AGROW>
        if anyAcceptedCount > 1, flags(end+1) = "orderAmbiguous"; end %#ok<AGROW>
        if numel(bestCand) >= 2 && bestCand(2).score >= bestCand(1).score - 0.2
            flags(end+1) = "ambiguous"; %#ok<AGROW>
        end
        bestName = ""; bestScore = NaN;
        if ~isempty(bestCand)
            bestName = string(bestCand(1).name); bestScore = bestCand(1).score;
        end
        rows(end+1, :) = {bestH.c, s.cSigma * bestH.scale, bestH.orders, ...
            ttF(s.members), bestName, bestScore, {bestCand}, {flags}}; %#ok<AGROW>
    end

    % Singles: try every order against the database; best (l, candidate) wins.
    for j = singleIdx(:).'
        bestCand = []; bestC = NaN; bestL = NaN;
        for l = 1:opts.MaxOrder
            cand = rankCandidates(films, dF(j) * l, sub.a, opts.WindowPad);
            if ~isempty(cand) && (isempty(bestCand) || cand(1).score > bestCand(1).score)
                bestCand = cand; bestC = dF(j) * l; bestL = l;
            end
        end
        if isempty(bestCand)
            unassigned(end+1) = ttF(j); %#ok<AGROW>
        else
            flags = "singlePeak";
            if numel(bestCand) >= 2 && bestCand(2).score >= bestCand(1).score - 0.2
                flags(end+1) = "ambiguous"; %#ok<AGROW>
            end
            rows(end+1, :) = {bestC, NaN, bestL, ttF(j), ...
                string(bestCand(1).name), bestCand(1).score, {bestCand}, {flags}}; %#ok<AGROW>
        end
    end

    R.series = cell2table(rows, 'VariableNames', ...
        {'cMeas','cSigma','orders','twoTheta','bestMatch','bestScore', ...
         'candidates','flags'});
    R.unassigned = unassigned(:);
    R.ghosts = ghosts;
    R.lambda = lambda;
    R.notes = strings(0, 1);
    if any(cellfun(@(c) any(string({c.name}) == "PZT"), R.series.candidates))
        R.notes(end+1) = PZT_CAVEAT;
    end
end

% ====================== local functions ======================

function [tt, counts] = normalizeInput(peaks)
    if isnumeric(peaks)
        tt = peaks(:); counts = nan(size(tt));
    elseif isstruct(peaks) && numel(peaks) > 1 || ...
           (isstruct(peaks) && isscalar(peaks) && isscalar(peaks.twoTheta))
        tt = [peaks.twoTheta].'; counts = [peaks.counts].';
    elseif isstruct(peaks) && isscalar(peaks)   % scan struct
        pk = xrdc.peaks.findPeaks(peaks, ...
            'MinProminence', max(peaks.counts) * 0.05);
        if isempty(pk), tt = []; counts = []; return, end
        tt = [pk.twoTheta].'; counts = [pk.counts].';
    else
        tt = []; counts = [];
    end
end

function m = meanOrNaN(v)
    if isempty(v), m = NaN; else, m = mean(v); end
end

function hyps = orderHypotheses(s, maxOrder)
%ORDERHYPOTHESES  The series as grouped, doubled, and (if even-only) halved.
    hyps = struct('kind', "asIs", 'c', s.c, 'orders', s.orders, 'scale', 1);
    if 2 * max(s.orders) <= maxOrder
        hyps(end+1) = struct('kind', "orderDoubled", 'c', 2 * s.c, ...
            'orders', 2 * s.orders, 'scale', 2);
    end
    if s.evenOnly
        hyps(end+1) = struct('kind', "orderHalved", 'c', s.c / 2, ...
            'orders', s.orders / 2, 'scale', 0.5);
    end
end

function cand = rankCandidates(films, cMeas, aSub, pad)
%RANKCANDIDATES  Accepted candidates for a measured c, best first.
    cand = struct('name', {}, 'score', {}, 'misfit', {}, 'cBulk', {}, ...
        'cPred', {}, 'strainVsBulk', {}, 'relaxation', {}, 'x', {}, ...
        'hasComposition', {});
    for e = films
        [hull, cBulk, cPred, x] = candidateWindow(e, cMeas, aSub);
        if cMeas >= hull(1) && cMeas <= hull(2)
            misfit = 0;
        else
            misfit = min(abs(cMeas - hull)) / e.c;
        end
        if misfit > pad, continue, end
        denom = cBulk - cPred;
        if abs(denom) < 1e-4, relax = NaN;
        else, relax = (cMeas - cPred) / denom; end
        cand(end+1) = struct('name', string(e.name), ...
            'score', max(0, 1 - misfit / pad), 'misfit', misfit, ...
            'cBulk', cBulk, 'cPred', cPred, ...
            'strainVsBulk', (cMeas - cBulk) / cBulk, ...
            'relaxation', relax, 'x', x, ...
            'hasComposition', ~isempty(e.composition)); %#ok<AGROW>
    end
    if isempty(cand), return, end
    % Rank: score desc, |cMeas-cPred| asc, fixed-composition first (parsimony).
    key = [-[cand.score]; abs(cMeas - [cand.cPred]); [cand.hasComposition]].';
    [~, order] = sortrows(key);
    cand = cand(order);
end

function [hull, cBulk, cPred, x] = candidateWindow(e, cMeas, aSub)
%CANDIDATEWINDOW  [min,max] of {c_bulk, c_pred} (over composition for PZT).
    f = strainFactor(e);
    x = NaN;
    if isempty(e.composition)
        cBulk = e.c;
        epsPar = (aSub - e.a) / e.a;
        cPred  = e.c * (1 - f * epsPar);
        hull   = [min(cBulk, cPred), max(cBulk, cPred)];
    else
        xg = linspace(min(e.composition.x), max(e.composition.x), 105);
        ab = interp1(e.composition.x, e.composition.a, xg);
        cb = interp1(e.composition.x, e.composition.c, xg);
        cp = cb .* (1 - f * (aSub - ab) ./ ab);
        hull = [min([cb cp]), max([cb cp])];
        % pseudomorphic inversion: cp is monotonic in x (a(x) increasing)
        if cMeas >= min(cp) && cMeas <= max(cp)
            x = interp1(cp, xg, cMeas);
            cBulk = interp1(xg, cb, x);
            cPred = cMeas;                       % on the pseudomorphic line
        else
            cBulk = e.c;                          % nominal entry
            epsPar = (aSub - e.a) / e.a;
            cPred  = e.c * (1 - f * epsPar);
        end
    end
end

function f = strainFactor(e)
    if isfield(e.elastic, 'nu') && ~isempty(e.elastic.nu)
        f = 2 * e.elastic.nu / (1 - e.elastic.nu);
    else
        f = 2 * e.elastic.c13 / e.elastic.c33;
    end
end
```

Implementation notes for the engineer:
- `normalizeInput` distinguishes a findPeaks struct ARRAY (scalar `.twoTheta` per element) from a scan struct (vector `.twoTheta`). The condition shown handles a 1-element findPeaks array via the `isscalar(peaks.twoTheta)` check.
- In `candidateWindow`'s pseudomorphic-inversion branch, `cPred = cMeas` makes the relaxation come out 0 by construction for an on-the-line film — intended: under the pseudomorphic assumption that point IS the prediction at the inferred x.
- Struct-array literal assignments (`hyps(end+1) = struct(...)`) require identical field orders — keep them exactly as written.
- If `sortrows` on the mixed key array misbehaves with logicals, cast: `double([cand.hasComposition])`.

- [ ] **Step 4: Run tests, fix until green**

`matlab -batch "results = runtests('tests/testIdentify.m'); disp(table(results)); assertSuccess(results)"`
Expected: all PASS. The synthetic numbers in the tests were derived analytically; if a tolerance assertion fails by a hair, re-derive the expected value rather than loosening the tolerance.

- [ ] **Step 5: Run the FULL suite** (regression gate)

From `xrdc-matlab/`: `matlab -batch "results = runtests; disp(table(results)); assertSuccess(results)"`
Expected: everything green (skips for gated real-data/toolbox tests are fine).

- [ ] **Step 6: Commit**

```bash
git add xrdc-matlab/+xrdc/+lattice/identifyMaterial.m xrdc-matlab/tests/testIdentify.m
git commit -m "feat: strain-aware material identification from (00l) peaks"
```

---

### Task 5: Real-data validation tests (S25 + S31)

**Files:**
- Modify: `xrdc-matlab/tests/testIdentify.m` (append)

Two gated tests against Tushar's scans (paths relative to repo root; the tests run from `xrdc-matlab/`, hence the `..`). S25 = SRO film on STO; S31 = PTO on SRO on STO (three phases in one scan).

- [ ] **Step 1: Write the tests** (append to `tests/testIdentify.m`)

```matlab
% ---------- real-data validation (gated) ----------

function testIdentifyS25RealData(tc)
    p = fullfile('..', 'validation', 'tushar', 'input', ...
        'TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_2theta_omega_05062026.txt');
    tc.assumeTrue(isfile(p), 'S25 validation scan not present');
    scan = xrdc.io.readScan(p);
    R = xrdc.lattice.identifyMaterial(scan, 1.5406, Substrate="SrTiO3");
    tc.verifyTrue(R.substrate.found);
    tc.verifyTrue(any(R.series.bestMatch == "SrRuO3"), ...
        'SRO film series not identified');
end

function testIdentifyS31Heterostructure(tc)
    p = fullfile('..', 'validation', 'tushar', 'input', ...
        'Heterostructure raw data', ...
        ['TR_S31_1_PTO on SRO_STO(100)_580c_150mT_and_700c_ 200mT', ...
         '_10500sh_5hz_2theta_omega_05202026.txt']);
    tc.assumeTrue(isfile(p), 'S31 validation scan not present');
    scan = xrdc.io.readScan(p);
    R = xrdc.lattice.identifyMaterial(scan, 1.5406, Substrate="SrTiO3");
    tc.verifyTrue(R.substrate.found);
    % PTO must appear among candidates of some series (bestMatch may be
    % flagged ambiguous vs dilute PZT - that is correct behaviour).
    hasPTO = any(cellfun(@(c) ~isempty(c) && ...
        any(string({c.name}) == "PbTiO3"), R.series.candidates));
    tc.verifyTrue(hasPTO, 'PTO not among any series candidates');
end
```

CAUTION: the S31 filename contains a space before `200mT` (`700c_ 200mT`) — that is faithful to the file on disk; verify with `dir` before assuming a typo. If `readScan`'s default peak prominence misses the film peaks (films are weak next to substrate), pass explicit peaks instead: `pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts)*0.01); R = xrdc.lattice.identifyMaterial(pk, ...)` — adjust the test accordingly and note it in a comment.

- [ ] **Step 2: Run the tests** — same testIdentify command. Expected: PASS (these run for real on this machine; the data is present).
- [ ] **Step 3: If S25/S31 fail, debug for real.** This is the validation moment of the whole feature; do not weaken assertions to pass. Inspect with `R = xrdc.lattice.identifyMaterial(...)`, look at `R.series.cMeas` vs expectations (SRO film c ≈ 3.95 Å; PTO film c ≈ 4.15 Å), and tune ONLY: findPeaks prominence in the test, or document a justified default change.
- [ ] **Step 4: Commit**

```bash
git add xrdc-matlab/tests/testIdentify.m
git commit -m "test: validate material ID on S25 (SRO/STO) and S31 (PTO/SRO/STO)"
```

---

### Task 6: Standalone build + docs

**Files:**
- Modify: `xrdc-matlab/build/buildStandalone.m:96`
- Modify: `xrdc-matlab/docs/FEATURES.md`

- [ ] **Step 1: Bundle the database.** In `buildStandalone.m` change line 96:

```matlab
    for f = ["substrates.json", "xrayLines.json", "materials.json"]
```

- [ ] **Step 2: Document the feature.** In `docs/FEATURES.md`, add under the θ-2θ capability section (near the line `- ✅ θ-2θ (twoThetaOmega) — peaks + substrate overlay`):

```markdown
- ✅ Material identification from (00l) peaks — `xrdc.lattice.identifyMaterial`: harmonic-series grouping, Kβ/W-Lα ghost filter, strain-aware candidate matching (pseudomorphic c prediction), ranked candidate sets, PZT composition estimate (strain-caveated). Database: `+xrdc/+data/materials.json` (STO, SRO, PTO, PZT). **[validated vs S25 + S31]**
```

- [ ] **Step 3: Run full suite** (`matlab -batch "results = runtests; assertSuccess(results)"` from `xrdc-matlab/`). Expected: green.
- [ ] **Step 4: Commit**

```bash
git add xrdc-matlab/build/buildStandalone.m xrdc-matlab/docs/FEATURES.md
git commit -m "build+docs: bundle materials.json, document material ID"
```

---

### Task 7: GUI wiring

**Files:**
- Modify: `xrdc-matlab/xrdcApp.m` — `buildAnalysisPanel` (case `'twothetaomega'`, ~line 311), new helpers after `addCheck` (~line 365), new callback + report/annotation functions near `runThetaTwoTheta` (~line 618)

Presentation note: the spec's "results dialog" is implemented as the existing Results panel (the monospace `resultsArea` textarea every other analysis uses) plus on-plot labels — consistent with the app; nested ranked candidates don't fit a `uitable`. This divergence is recorded in the spec amendment in Step 6.

- [ ] **Step 1: Add controls to the θ-2θ panel.** In `buildAnalysisPanel`, replace the `case 'twothetaomega'` body:

```matlab
        case 'twothetaomega'
            row = addEdit (g, row, 'Min prom (%)',   '5',   @(v) onParamChange(fig, 'promPct',   v));
            subs = identifiableSubstrates();
            row = addDrop (g, row, 'Substrate', subs, 'SrTiO3', ...
                                                      @(v) onParamChange(fig, 'substrate', v));
            row = addIdentifyButton(g, row, fig); %#ok<NASGU>
```

- [ ] **Step 2: Add the helpers** (after `addCheck`, before `appTheme`):

```matlab
function subs = identifiableSubstrates()
%IDENTIFIABLESUBSTRATES  Dropdown items: materials usable as a substrate.
    M = xrdc.lattice.loadMaterials();
    ok = ismember(string({M.role}), ["substrate", "both"]);
    subs = cellstr(string({M(ok).name}));
end

function row = addIdentifyButton(g, row, fig)
    T = appTheme();
    b = uibutton(g, 'Text', 'Identify Material', ...
        'FontName', T.font, 'FontSize', 12, ...
        'BackgroundColor', T.btn, 'FontColor', T.text, ...
        'ButtonPushedFcn', @(~,~) onIdentifyMaterial(fig));
    b.Layout.Row = row; b.Layout.Column = [1 2];
    row = row + 1;
end
```

- [ ] **Step 3: Add the callback + report + annotation** (place after `runThetaTwoTheta`):

```matlab
function onIdentifyMaterial(fig)
%ONIDENTIFYMATERIAL  Run material ID on the current theta-2theta peaks.
    st = fig.UserData;
    if isempty(st.scan), return, end
    promPct = getNum(st.params, 'promPct', 5);
    pk = xrdc.peaks.findPeaks(st.scan, ...
        'MinProminence', max(st.scan.counts) * promPct / 100);
    if isempty(pk)
        uialert(fig, 'No peaks detected - lower "Min prom (%)" and retry.', ...
            'Identify material', 'Icon', 'warning');
        return
    end
    sub = getStr(st.params, 'substrate', 'SrTiO3');
    lambda = 1.5406;
    if isfield(st.scan, 'lambda') && ~isempty(st.scan.lambda) ...
            && isfinite(st.scan.lambda)
        lambda = st.scan.lambda;
    end
    try
        R = xrdc.lattice.identifyMaterial(pk, lambda, Substrate=string(sub));
    catch ME
        uialert(fig, sprintf('Identification failed:\n\n%s', ME.message), ...
            'Identify material', 'Icon', 'error');
        return
    end
    annotateIdentification(st.ax, R);
    writeResults(fig, identificationReport(R));
end

function annotateIdentification(ax, R)
%ANNOTATEIDENTIFICATION  Material + (00l) labels above identified peaks.
    hold(ax, 'on');
    yl = ylim(ax);
    for i = 1:numel(R.substrate.twoTheta)
        labelPeak(ax, R.substrate.twoTheta(i), yl, ...
            sprintf('%s', shortName(R.substrate.name)), [0 0 0]);
    end
    for i = 1:height(R.series)
        name = R.series.bestMatch(i);
        if name == "", name = "?"; end
        tts = R.series.twoTheta{i};  ords = R.series.orders{i};
        for j = 1:numel(tts)
            labelPeak(ax, tts(j), yl, ...
                sprintf('%s (00%d)', shortName(name), ords(j)), ...
                [0.65 0.15 0.15]);
        end
    end
    hold(ax, 'off');
end

function labelPeak(ax, tt, yl, txt, color)
    text(ax, tt, yl(2) * 0.7, txt, 'Rotation', 90, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'Color', color, 'Interpreter', 'none');
end

function s = shortName(name)
%SHORTNAME  Compact display alias for plot labels.
    map = struct('SrTiO3', "STO", 'SrRuO3', "SRO", 'PbTiO3', "PTO");
    n = char(name);
    if isfield(map, n), s = map.(n); else, s = string(name); end
end

function lines = identificationReport(R)
%IDENTIFICATIONREPORT  Results-panel text for an identifyMaterial run.
    lines = {};
    if R.substrate.found
        lines{end+1} = sprintf('Substrate %s: confirmed (%d peaks), c = %.4f A', ...
            R.substrate.name, numel(R.substrate.twoTheta), R.substrate.cMeas);
    else
        lines{end+1} = sprintf('Substrate %s: NOT FOUND in scan', R.substrate.name);
    end
    if ~isempty(R.ghosts) && height(R.ghosts) > 0
        lines{end+1} = sprintf('Filtered %d Kbeta/W-La ghost peak(s)', height(R.ghosts));
    end
    for i = 1:height(R.series)
        lines{end+1} = '';
        lines{end+1} = sprintf('Series %d: c = %.4f A (orders: %s)', ...
            i, R.series.cMeas(i), strjoin(string(R.series.orders{i}), ','));
        cand = R.series.candidates{i};
        if isempty(cand)
            lines{end+1} = '  -> no database match (unidentified)';
        end
        for k = 1:numel(cand)
            extra = '';
            if ~isnan(cand(k).x)
                extra = sprintf('   x(Zr) ~ %.2f', cand(k).x);
            end
            lines{end+1} = sprintf( ...
                '  -> %-7s score %.2f   strain vs bulk %+.2f%%   relax %.2f%s', ...
                cand(k).name, cand(k).score, 100 * cand(k).strainVsBulk, ...
                cand(k).relaxation, extra);
        end
        fl = R.series.flags{i};
        if ~isempty(fl)
            lines{end+1} = sprintf('  flags: %s', strjoin(fl, ', '));
        end
    end
    if ~isempty(R.unassigned)
        lines{end+1} = '';
        lines{end+1} = sprintf('Unassigned peaks: %s', ...
            strjoin(compose('%.2f', R.unassigned), ', '));
    end
    for nt = R.notes(:).'
        lines{end+1} = '';
        lines{end+1} = char("NOTE: " + nt);
    end
end
```

- [ ] **Step 4: Run the full suite** — the GUI has no automated tests, but the suite catches syntax errors in `xrdcApp.m`? It does NOT (xrdcApp is not under test). Instead lint it: `matlab -batch "issues = codeIssues('xrdcApp.m'); disp(issues.Issues); assert(~any(issues.Issues.Severity == 'error'))"` from `xrdc-matlab/`. Expected: no errors (info/warnings acceptable if pre-existing).
- [ ] **Step 5: Manual smoke test** (requires display): `matlab -r xrdcApp`, load `..\validation\tushar\input\TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_2theta_omega_05062026.txt`, click **Identify Material**. Verify: STO labels on substrate peaks, SRO on film peaks, ranked candidates + flags in the Results panel. Screenshot for the session notes if convenient.
- [ ] **Step 6: Amend the spec presentation line.** In `docs/superpowers/specs/2026-06-09-material-id-design.md`, GUI section: replace "shows the ranked results table in a dialog" with "shows the ranked results in the Analysis panel's Results area (the app-wide pattern; nested ranked candidates don't fit a uitable)".
- [ ] **Step 7: Commit**

```bash
git add xrdc-matlab/xrdcApp.m docs/superpowers/specs/2026-06-09-material-id-design.md
git commit -m "feat(gui): Identify Material action for theta-2theta scans"
```

---

## Self-review notes (already applied)

- **Spec coverage:** database (T1), ghost filter (T2), grouping + c/2 guard (T3+T4 hypotheses), substrate-as-input + strain-aware naming + ranked candidates + relaxation + PZT caveat (T4), real-data + synthetic tests (T4+T5), build bundle + docs (T6), GUI (T7). The spec's "dialog" is amended in T7 Step 6.
- **Order-ambiguity direction:** the spec wrote the guard as c-vs-c/2; the greedy grouper actually lands on the SMALLEST consistent c, so the common correction is doubling (2c), not halving. T4's `orderHypotheses` evaluates asIs/doubled/halved and the test `testIdentifyOrderDoubling` pins the behaviour. This is a refinement of the spec's mechanism, same intent.
- **Type consistency:** `loadMaterials` lookup signature used by T4 and T7; `groupHarmonicSeries` outputs `.members/.orders/.c/.cSigma/.evenOnly` consumed in T4; `filterGhostPeaks` `[keep, ghosts-table]` consumed in T4; `R.series` columns consumed by T7's report/annotation. Checked.
