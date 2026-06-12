# Auto Min-Prominence Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** `xrdc.peaks.findPeaks` auto-selects a peak criterion that finds the physically correct peak set on θ-2θ scans (substrate + film, no noise/duplicates) without manual threshold tuning; the GUI and `identifyMaterial` use it by default.

**Architecture:** When `MinProminence` is omitted, detect on `log10(counts)` with a 0.3-decade prominence threshold (scale-invariant across the ~6-decade dynamic range), reject candidates that fail a 5σ Poisson test against the local `movmedian` background, then report linear-domain metrics from a second `findpeaks` pass intersected by index (log10 is monotone, so maxima indices coincide). Explicit numeric `MinProminence` keeps the existing linear path byte-for-byte.

**Tech Stack:** MATLAB R2022b+, `arguments` blocks, Signal Processing Toolbox `findpeaks` with the existing `xrdc.peaks.findpeaks_fallback` fallback (it already supports every option the auto path uses: `MinPeakProminence`, `MinPeakDistance`, `MinPeakHeight`, `WidthReference` — the last is parsed and unused there, which is fine because widths are only consumed from the *linear* pass where the fallback computes half-height widths). `movmedian` is base MATLAB.

**Spec:** `docs/superpowers/specs/2026-06-12-auto-prominence-design.md`

**House rules:** run the FULL suite (`matlab -batch "results = runtests; disp(table(results)); assertSuccess(results)"` from `xrdc-matlab/`) before every commit. Commit after each green step. NEVER push.

---

## File map

- Modify: `xrdc-matlab/+xrdc/+peaks/findPeaks.m` — auto path (Task 1)
- Modify: `xrdc-matlab/tests/testPeaks.m` — new tests (Tasks 1, 2)
- Modify: `xrdc-matlab/+xrdc/+lattice/identifyMaterial.m:184-188` — scan-struct input (Task 3)
- Modify: `xrdc-matlab/tests/testIdentify.m` — scan-struct test (Task 3)
- Modify: `xrdc-matlab/xrdcApp.m:312,639-675` — GUI auto default + shared helper (Task 4)
- Create: `validation/runAutoProminence.m` — calibration sweep script (Task 5)
- Modify: `xrdc-matlab/tests/testPeaks.m` — gated real-data tests (Task 5)
- Modify: `xrdc-matlab/docs/USER_GUIDE.md`, `xrdc-matlab/docs/SCIENTIFIC_ASSUMPTIONS.md`, `xrdc-matlab/docs/FEATURES.md`, `validation/DATA_SWEEP.md` (Task 6)

All MATLAB commands run from `xrdc-matlab/` unless stated otherwise.

---

### Task 1: Log-domain auto detection in `findPeaks`

**Files:**
- Modify: `xrdc-matlab/+xrdc/+peaks/findPeaks.m`
- Test: `xrdc-matlab/tests/testPeaks.m`

- [ ] **Step 1: Write three failing tests**

Append to `xrdc-matlab/tests/testPeaks.m`, directly after `testFindPeaksBadScan` (before the `findPeaksLegacy` section header):

```matlab
% ---------- findPeaks auto prominence (log-domain) ----------

function testFindPeaksAutoWideDynamicRange(tc)
    % Substrate-like 1e6 peak + film-like 1.5e3 peak on a 30-count
    % background. No single linear threshold separates both from noise;
    % the log-domain auto criterion must find exactly the two.
    x = (20:0.01:50).';
    scan = syntheticScan(x, [22.7, 46.5], [0.05, 0.30], ...
                              [1e6, 1.5e3], "gauss", 30);
    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyLength(pk, 2);
    tc.verifyEqual(sort([pk.twoTheta]), [22.7, 46.5], 'AbsTol', 0.05);
end

function testFindPeaksAutoRejectsQuantisationNoise(tc)
    % Near-zero-count baseline: 0..3-count jitter has ~0.48 decades of
    % log10 "prominence" but is statistically nothing. The Poisson guard
    % must reject every candidate.
    x = (20:0.02:40).';
    y = double(mod(0:numel(x)-1, 4)).';   % 0 1 2 3 0 1 2 3 ...
    scan = xrdc.io.emptyScan();
    scan.twoTheta = x;
    scan.counts   = y;
    scan.sourceFormat = "synthetic";
    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyEmpty(pk);
end

function testFindPeaksAutoSuppressesFringeRipple(tc)
    % ±10% multiplicative ripple (thickness-fringe-like) is ~0.08 decades
    % in log space — below the 0.3-decade threshold everywhere — while the
    % main peak (~2 decades) survives. Exactly one peak.
    x = (25:0.01:35).';
    scan = syntheticScan(x, 30, 0.3, 1e4, "gauss", 100);
    scan.counts = scan.counts .* (1 + 0.10 * sin(2*pi*x/0.4));
    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyLength(pk, 1);
    tc.verifyEqual(pk(1).twoTheta, 30, 'AbsTol', 0.06);
end
```

- [ ] **Step 2: Run the new tests, verify all three FAIL**

Run: `matlab -batch "results = runtests('tests/testPeaks.m'); disp(table(results)); assertSuccess(results)"`

Expected failures under the current 1.5%-of-range default:
- `AutoWideDynamicRange`: threshold ≈ 15000 counts → film peak missed → 1 peak, not 2.
- `AutoRejectsQuantisationNoise`: threshold = max(1, 0.045) = 1 → every "3" is a peak → non-empty.
- `AutoSuppressesFringeRipple`: threshold ≈ 165 counts → ripple bumps on the peak flanks (linear prominence ~hundreds) pass → more than 1 peak.

If any of the three passes, STOP — the test is wrong; fix it before implementing.

- [ ] **Step 3: Implement the auto path**

In `xrdc-matlab/+xrdc/+peaks/findPeaks.m`:

**(a)** Replace the `'MinProminence'` docstring block (lines 17–18):

```matlab
%     'MinProminence'  (default: automatic, log-domain)
%         Required prominence above neighbouring troughs, in counts.
%         When omitted, an automatic scale-invariant criterion is used:
%         a peak must rise >= 0.3 decades above its neighbouring troughs
%         on the log10(counts) curve AND stand >= 5 Poisson sigma above
%         the local (~1 degree) median background. One default covers a
%         1e6-count substrate peak and a 1e3-count film peak in the same
%         scan. Counts below 1 are clamped before the log, so sub-1-cps
%         baselines yield no peaks. Pass a numeric value (counts) for the
%         classic fixed linear threshold.
```

**(b)** Replace lines 83–128 — everything from the `% Default prominence = 1.5% ...` comment through the `findpeaks`/fallback call (keep the window-cropping and `step`/`minDistSamples` code in place; the listing below shows the full region after the edit so surrounding context is unambiguous):

```matlab
    % Restrict to requested 2θ window before calling findpeaks so the
    % MinPeakDistance constraint is interpreted in-window.
    inWin = x >= options.TwoThetaRange(1) & x <= options.TwoThetaRange(2);
    if ~any(inWin)
        peaks = emptyPeakArray();
        return
    end
    xw = x(inWin);
    yw = y(inWin);
    idxInScan = find(inWin);

    % MinPeakDistance is in samples (findpeaks semantics); convert from
    % 2θ using the median step. For irregular spacing this is still a
    % reasonable default.
    step = median(diff(xw));
    if step <= 0 || ~isfinite(step)
        minDistSamples = 1;
    else
        minDistSamples = max(1, round(options.MinSeparation / step));
    end

    if isnan(options.MinProminence)
        % Automatic log-domain criterion (see help text above).
        [pks, locs, widths, proms] = autoDetect(yw, step, minDistSamples, options);
    else
        fpArgs = { ...
            'MinPeakProminence', options.MinProminence, ...
            'MinPeakHeight',     options.MinHeight, ...
            'MinPeakDistance',   minDistSamples, ...
            'WidthReference',    char(options.WidthReference)};
        if options.MinWidth > 0
            fpArgs = [fpArgs, {'MinPeakWidth', options.MinWidth / max(step,eps)}];
        end
        if isfinite(options.MaxWidth)
            fpArgs = [fpArgs, {'MaxPeakWidth', options.MaxWidth / max(step,eps)}];
        end
        [pks, locs, widths, proms] = callFindpeaks(yw, fpArgs{:});
    end
```

Note: the old `minProm = options.MinProminence; if isnan(minProm) ...` block (former lines 83–88) is deleted; the manual branch passes `options.MinProminence` directly.

**(c)** Add two local functions at the end of the file, before `blankPeak()`:

```matlab
% -------------------------------------------------------------------------

function [pks, locs, widths, proms] = autoDetect(yw, step, minDistSamples, options)
%AUTODETECT  Log-domain auto prominence (MinProminence omitted).
%   A peak must (a) rise PROM_DECADES above its neighbouring troughs on
%   the log10 intensity curve — scale-invariant across the ~6-decade
%   dynamic range of a θ-2θ scan — and (b) clear a Poisson significance
%   test against the local median background, which rejects the large
%   fake log-prominence of quantisation jitter on near-zero-count
%   baselines. Reported metrics stay linear-domain: log10 (clamped at 1)
%   is monotone non-decreasing, so every accepted log-domain maximum is a
%   linear local maximum at the same index.
    PROM_DECADES  = 0.3;   % ≈2x above the surrounding troughs
    NOISE_SIGMAS  = 5;     % Poisson significance vs local background
    BG_WINDOW_DEG = 1.0;   % local-background median window, in 2θ

    pks = []; locs = []; widths = []; proms = [];

    yLog = log10(max(yw, 1));
    [~, locsLog] = callFindpeaks(yLog, ...
        'MinPeakProminence', PROM_DECADES, ...
        'MinPeakDistance',   minDistSamples);
    if isempty(locsLog), return, end

    if step > 0 && isfinite(step)
        bgWin = max(3, round(BG_WINDOW_DEG / step));
    else
        bgWin = min(numel(yw), 51);
    end
    bg  = movmedian(yw, bgWin);
    sig = (yw(locsLog) - bg(locsLog)) >= NOISE_SIGMAS * sqrt(max(bg(locsLog), 1));
    locsKeep = locsLog(sig);
    if isempty(locsKeep), return, end

    % Linear metrics: enumerate all linear local maxima with widths and
    % prominences, then keep the accepted subset by index.
    [pksAll, locsAll, widthsAll, promsAll] = callFindpeaks(yw, ...
        'MinPeakHeight',  options.MinHeight, ...
        'WidthReference', char(options.WidthReference));
    [locs, ia] = intersect(locsAll, locsKeep);
    pks = pksAll(ia); widths = widthsAll(ia); proms = promsAll(ia);

    keep = true(size(locs));
    if options.MinWidth > 0
        keep = keep & widths >= options.MinWidth / max(step, eps);
    end
    if isfinite(options.MaxWidth)
        keep = keep & widths <= options.MaxWidth / max(step, eps);
    end
    pks = pks(keep); locs = locs(keep); widths = widths(keep); proms = proms(keep);
end

function [pks, locs, widths, proms] = callFindpeaks(y, varargin)
%CALLFINDPEAKS  Signal Processing Toolbox findpeaks, or the fallback.
    if isempty(which('findpeaks'))
        [pks, locs, widths, proms] = xrdc.peaks.findpeaks_fallback(y, varargin{:});
    else
        [pks, locs, widths, proms] = findpeaks(y, varargin{:});
    end
end
```

Also update the stale top-of-file reference: the existing comment block `% Use Signal Processing Toolbox findpeaks if available; otherwise fallback.` above the old call site is superseded by `callFindpeaks` — delete it if it survives the (b) replacement.

- [ ] **Step 4: Run testPeaks, verify the three new tests PASS and the pre-existing ones still pass**

Run: `matlab -batch "results = runtests('tests/testPeaks.m'); disp(table(results)); assertSuccess(results)"`
Expected: PASS, including `testFindPeaksSingleLorentz` / `MultiplePeaks` / `HeightThreshold` / `RangeCrop` (these exercise the default path and now go through `autoDetect`).

- [ ] **Step 5: Run the FULL suite**

Run: `matlab -batch "results = runtests; disp(table(results)); assertSuccess(results)"`
Expected: green. `testXrr` / `testIdentify` / examples-driven tests pass explicit prominences and must be unaffected.

- [ ] **Step 6: Commit**

```bash
git add xrdc-matlab/+xrdc/+peaks/findPeaks.m xrdc-matlab/tests/testPeaks.m
git commit -m "feat(peaks): log-domain auto prominence with Poisson noise guard"
```

---

### Task 2: Option-interplay and backward-compat tests

**Files:**
- Test: `xrdc-matlab/tests/testPeaks.m`

These are regression specs; some may pass immediately (that is the point — they pin behaviour). Any failure means Task 1's implementation is wrong: fix `findPeaks.m`, not the test.

- [ ] **Step 1: Append three tests** (after `testFindPeaksAutoSuppressesFringeRipple`):

```matlab
function testFindPeaksExplicitProminenceUnchanged(tc)
    % Explicit MinProminence keeps the classic fixed linear path: a
    % 100-count bump fails a 300-count threshold even though the auto
    % criterion would accept it.
    x = (20:0.02:40).';
    scan = syntheticScan(x, [25, 35], [0.3, 0.3], [100, 5000], ...
                              "lorentz", 50);
    pk = xrdc.peaks.findPeaks(scan, 'MinProminence', 300);
    tc.verifyLength(pk, 1);
    tc.verifyEqual(pk(1).twoTheta, 35, 'AbsTol', 0.05);
end

function testFindPeaksAutoRespectsMinHeight(tc)
    % MinHeight applies to the linear pass in auto mode.
    x = (20:0.02:40).';
    scan = syntheticScan(x, [25, 35], [0.3, 0.3], [500, 5000], ...
                              "lorentz", 50);
    pk = xrdc.peaks.findPeaks(scan, 'MinHeight', 2000);
    tc.verifyLength(pk, 1);
    tc.verifyEqual(pk(1).twoTheta, 35, 'AbsTol', 0.05);
end

function testFindPeaksAutoRespectsMinSeparation(tc)
    % Kα-split-like pair 0.12° apart merges to the taller member when
    % MinSeparation exceeds the split (auto mode).
    x = (40:0.005:50).';
    scan = syntheticScan(x, [46.50, 46.62], [0.05, 0.05], ...
                              [1e6, 5e5], "gauss", 100);
    pk = xrdc.peaks.findPeaks(scan, 'MinSeparation', 0.2);
    tc.verifyLength(pk, 1);
    tc.verifyEqual(pk(1).twoTheta, 46.50, 'AbsTol', 0.02);
end
```

- [ ] **Step 2: Run testPeaks**

Run: `matlab -batch "results = runtests('tests/testPeaks.m'); disp(table(results)); assertSuccess(results)"`
Expected: PASS. On failure, diagnose against the auto-path semantics in Task 1 Step 3 and fix the implementation.

- [ ] **Step 3: Run the FULL suite** — same command as Task 1 Step 5. Expected: green.

- [ ] **Step 4: Commit**

```bash
git add xrdc-matlab/tests/testPeaks.m
git commit -m "test(peaks): pin auto-mode option interplay and explicit-prominence path"
```

---

### Task 3: `identifyMaterial` scan-struct input uses auto detection

**Files:**
- Modify: `xrdc-matlab/+xrdc/+lattice/identifyMaterial.m:184-188` (`normalizeInput`), docstring line ~18
- Test: `xrdc-matlab/tests/testIdentify.m`

- [ ] **Step 1: Write the failing test**

Append to `xrdc-matlab/tests/testIdentify.m`, just above the `% ---------- real-data validation (gated) ----------` header:

```matlab
function testIdentifyScanStructAutoFindsWeakFilm(tc)
    % Scan-struct input runs auto peak detection. A film 3 decades below
    % the substrate (0.1% of max counts — far under the old hardcoded 5%
    % rule) must still be detected and identified.
    lambda = 1.5406;
    x = (20:0.01:75).';
    y = 30 * ones(size(x));
    cSTO = 3.905; cPTO = 4.1511;
    for l = 1:3
        ttS = 2*asind(l*lambda/(2*cSTO));
        ttP = 2*asind(l*lambda/(2*cPTO));
        y = y + 2e6 * exp(-(x - ttS).^2 / (2*0.02^2));
        y = y + 2e3 * exp(-(x - ttP).^2 / (2*0.05^2));
    end
    scan = xrdc.io.emptyScan();
    scan.twoTheta = x;
    scan.counts   = y;
    R = xrdc.lattice.identifyMaterial(scan, lambda, Substrate="SrTiO3");
    tc.verifyTrue(R.substrate.found);
    tc.verifyTrue(any(R.series.bestMatch == "PbTiO3"), ...
        'weak PTO film series not identified from scan-struct input');
end
```

- [ ] **Step 2: Run testIdentify, verify the new test FAILS**

Run: `matlab -batch "results = runtests('tests/testIdentify.m'); disp(table(results)); assertSuccess(results)"`
Expected: FAIL — old 5% rule (threshold 1e5 counts) drops all PTO peaks, so no `PbTiO3` series.

- [ ] **Step 3: Change `normalizeInput`**

In `xrdc-matlab/+xrdc/+lattice/identifyMaterial.m`, replace:

```matlab
    elseif isstruct(peaks) && isscalar(peaks)   % scan struct
        pk = xrdc.peaks.findPeaks(peaks, ...
            'MinProminence', max(peaks.counts) * 0.05);
```

with:

```matlab
    elseif isstruct(peaks) && isscalar(peaks)   % scan struct
        % Auto log-domain prominence; 0.2° separation merges the Kα1/Kα2
        % substrate split before harmonic grouping.
        pk = xrdc.peaks.findPeaks(peaks, 'MinSeparation', 0.2);
```

And update the docstring input description (line ~18) from:

```matlab
%              - scan struct (vector .twoTheta/.counts) -> findPeaks is run
```

to:

```matlab
%              - scan struct (vector .twoTheta/.counts) -> findPeaks is
%                run with auto prominence and 0.2-deg separation
```

- [ ] **Step 4: Run testIdentify, verify PASS** — same command as Step 2. Expected: PASS, including all pre-existing synthetic tests (they pass 2θ vectors or peak structs, not scan structs, so they bypass `normalizeInput`'s scan branch).

- [ ] **Step 5: Run the FULL suite** — same command as Task 1 Step 5. Expected: green.

- [ ] **Step 6: Commit**

```bash
git add xrdc-matlab/+xrdc/+lattice/identifyMaterial.m xrdc-matlab/tests/testIdentify.m
git commit -m "feat(lattice): identifyMaterial scan-struct input uses auto peak detection"
```

---

### Task 4: GUI — `auto` default and shared θ-2θ peak helper

**Files:**
- Modify: `xrdc-matlab/xrdcApp.m` (three edits + one new local function)

The GUI has no automated test harness (established in the material-id plan); verification is lint + the full suite (catches nothing GUI-specific but guards the repo), then empirical checks in code review.

- [ ] **Step 1: Default the field to auto**

At `xrdcApp.m:312`, change:

```matlab
            row = addEdit (g, row, 'Min prom (%)',   '5',   @(v) onParamChange(fig, 'promPct',   v));
```

to:

```matlab
            row = addEdit (g, row, 'Min prom (%)',   'auto', @(v) onParamChange(fig, 'promPct',  v));
```

- [ ] **Step 2: Add the shared helper**

Insert directly above `function runThetaTwoTheta(fig)` (line 639):

```matlab
function pk = thetaPeaks(st)
%THETAPEAKS  Shared peak set for θ-2θ analyses (plot + Identify Material).
%   'Min prom (%)' numeric → manual threshold as % of max counts (the
%   classic behaviour). 'auto' (the default) or anything non-numeric →
%   findPeaks' log-domain auto criterion, with MinSeparation 0.2° so the
%   Kα1/Kα2 substrate split reports as one peak.
    raw = getStr(st.params, 'promPct', 'auto');
    pct = str2double(raw);
    if isfinite(pct)
        pk = xrdc.peaks.findPeaks(st.scan, ...
            'MinProminence', max(st.scan.counts) * pct / 100);
    else
        pk = xrdc.peaks.findPeaks(st.scan, 'MinSeparation', 0.2);
    end
end
```

- [ ] **Step 3: Use it in both call sites**

In `runThetaTwoTheta` replace:

```matlab
    promPct = getNum(st.params, 'promPct', 5);
    pk = xrdc.peaks.findPeaks(scan, ...
        'MinProminence', max(scan.counts) * promPct / 100);
```

with:

```matlab
    pk = thetaPeaks(st);
```

In `onIdentifyMaterial` replace:

```matlab
    promPct = getNum(st.params, 'promPct', 5);
    pk = xrdc.peaks.findPeaks(st.scan, ...
        'MinProminence', max(st.scan.counts) * promPct / 100);
    if isempty(pk)
        uialert(fig, 'No peaks detected - lower "Min prom (%)" and retry.', ...
            'Identify material', 'Icon', 'warning');
        return
    end
```

with:

```matlab
    pk = thetaPeaks(st);
    if isempty(pk)
        uialert(fig, ['No peaks detected - set "Min prom (%)" to auto, ' ...
            'or lower the percentage, and retry.'], ...
            'Identify material', 'Icon', 'warning');
        return
    end
```

(`runThetaTwoTheta` keeps its local `scan = st.scan` for plotting; only the peak call changes. The rocking-curve and XRR paths keep their explicit prominences — do not touch them.)

- [ ] **Step 4: Lint**

Run: `matlab -batch "issues = codeIssues('xrdcApp.m'); disp(issues.Issues); assert(~any(issues.Issues.Severity == 'error'))"`
Expected: no errors (pre-existing infos/warnings acceptable). Also confirm `getNum` is still referenced elsewhere in the file (it is — XRR/RSM params) so no orphan-helper warning appears.

- [ ] **Step 5: Run the FULL suite** — same command as Task 1 Step 5. Expected: green.

- [ ] **Step 6: Commit**

```bash
git add xrdc-matlab/xrdcApp.m
git commit -m "feat(gui): theta-2theta peak detection defaults to auto prominence"
```

---

### Task 5: Real-data calibration + gated regression tests

**Files:**
- Create: `validation/runAutoProminence.m`
- Test: `xrdc-matlab/tests/testPeaks.m`

Ground truth (from the design spec): substrate (00l) + film (00l) peaks only. STO 001/002/003 at 22.75/46.47/72.57°. A resolved Kα₂ companion of the STO 003 (~0.21° above it, beyond the 0.2° merge) is acceptable. No quantisation-noise or fringe picks; no over-segmentation.

- [ ] **Step 1: Write the calibration sweep script**

Create `validation/runAutoProminence.m`:

```matlab
%RUNAUTOPROMINENCE  Auto-prominence sweep over the theta-2theta data dump.
%   Prints, per scan: peak count, then 2-theta / counts / FWHM per peak.
%   Eyeball criteria (design spec 2026-06-12): STO scans show the 3
%   substrate peaks (22.75/46.47/72.57 +/- 0.05) plus film 00l peaks;
%   no two peaks closer than ~0.15 deg; no baseline-noise picks.
%   Run from the repo root with xrdc-matlab on the path:
%     cd xrdc-matlab && matlab -batch "addpath(pwd); cd ../validation; runAutoProminence"

names = { ...
    'TR_S04_PTO_STO(100)_750c_200mT_1000sh_3hz_2theta omega_04072026.txt', ...
    'TR_S05_PTO_STO(100)_600c_200mT_1000sh_2hz_2theta omega_04092026.txt', ...
    'TR_S06_PTO_LAO(100)_600c_200mT_1000sh_2hz_2 theta omega_04082026.txt', ...
    'TR_S07_PTO_STO(100)_600c_variable pressure_2000sh_3hz_2theta omega_04102026.txt', ...
    'TR_S08_PTO_STO(100)_550c_200mT_5000sh_3hz_2theta omega_04132026.txt', ...
    'TR_S10_PTO_STO(100)_500c_150mT_20000sh_5hz_2theta omega_04162026.txt', ...
    'TR_S11_PTO_STO(100)_580c_150mT_20000sh_5hz_2theta omega_04162026.txt'};

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
for k = 1:numel(names)
    p = fullfile(dataDir, names{k});
    if ~isfile(p)
        fprintf('\nSKIP (missing): %s\n', names{k});
        continue
    end
    scan = xrdc.io.readScan(p);
    pk = xrdc.peaks.findPeaks(scan, 'MinSeparation', 0.2);
    fprintf('\n%s\n  -> %d peak(s)\n', names{k}, numel(pk));
    for i = 1:numel(pk)
        fprintf('   %8.3f deg   %10.0f cts   fwhm %.3f deg\n', ...
            pk(i).twoTheta, pk(i).counts, pk(i).fwhm);
    end
end
```

- [ ] **Step 2: Run the sweep and judge every file against the eyeball criteria**

Run from `xrdc-matlab/`: `matlab -batch "addpath(pwd); cd ../validation; runAutoProminence"`

Record the full output in your task report. For each STO file (S04, S05, S07, S08, S10, S11): exactly one peak within 0.05° of each of 22.75/46.47/72.57; remaining peaks must be plausible film (00l) reflections (PTO c between bulk 4.15 Å and pseudomorphic ~4.26 Å puts 001/002/003 in roughly 20.8–21.4 / 42.5–43.6 / 66.0–67.7°) or the resolved STO-003 Kα₂; nothing in the low-count baseline. For S06 (LAO substrate): three substrate peaks near 23.4/48.0/74.2° plus film peaks, total well under 15.

- [ ] **Step 3: Recalibrate constants ONLY if Step 2 fails**

If a file shows noise picks or missing film peaks, adjust the constants in `autoDetect` (`findPeaks.m`) in this order of preference: `PROM_DECADES` within [0.25, 0.5], `NOISE_SIGMAS` within [4, 8], `BG_WINDOW_DEG` within [0.5, 2.0]. After ANY change: re-run `matlab -batch "results = runtests('tests/testPeaks.m'); disp(table(results)); assertSuccess(results)"` AND repeat Step 2 until both are clean. Record final constants and the reasoning in the task report. If no values in those ranges satisfy both, STOP and report BLOCKED with the conflicting evidence.

- [ ] **Step 4: Add gated real-data tests**

Append to `xrdc-matlab/tests/testPeaks.m` (end of file, before the `small helpers` section):

```matlab
% ---------- auto prominence: real data (gated on the local data dump) ----------

function d = dataDumpDir()
%Helper (not a test): repo-root data dump, anchored on this file's location
%   (the unittest runner cd's into tests/, so cwd-relative paths are fragile).
    d = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end

function pk = autoPeaks(tc, fname)
%Helper (not a test): gated read + auto peak detection, GUI-equivalent call.
    p = fullfile(dataDumpDir(), fname);
    tc.assumeTrue(isfile(p), sprintf('data dump scan not present: %s', fname));
    scan = xrdc.io.readScan(p);
    pk = xrdc.peaks.findPeaks(scan, 'MinSeparation', 0.2);
end

function verifySTOSubstrate(tc, pk)
%Helper (not a test): the STO 001/002/003 triplet must be present, once each.
    tt = [pk.twoTheta];
    for ref = [22.75, 46.47, 72.57]
        tc.verifyEqual(nnz(abs(tt - ref) < 0.05), 1, ...
            sprintf('expected exactly one peak within 0.05 deg of %.2f', ref));
    end
end

function testAutoProminenceS05RealScan(tc)
    % S05: previously needed a manually lowered threshold to see the film.
    pk = autoPeaks(tc, ...
        'TR_S05_PTO_STO(100)_600c_200mT_1000sh_2hz_2theta omega_04092026.txt');
    verifySTOSubstrate(tc, pk);
    tc.verifyGreaterThanOrEqual(numel(pk), 4, 'film peaks missing');
    tc.verifyLessThanOrEqual(numel(pk), 12, 'over-segmented');
    tc.verifyTrue(all(diff(sort([pk.twoTheta])) >= 0.15), ...
        'duplicate/split picks closer than 0.15 deg');
end

function testAutoProminenceS04RealScan(tc)
    pk = autoPeaks(tc, ...
        'TR_S04_PTO_STO(100)_750c_200mT_1000sh_3hz_2theta omega_04072026.txt');
    verifySTOSubstrate(tc, pk);
    tc.verifyGreaterThanOrEqual(numel(pk), 4, 'film peaks missing');
    tc.verifyLessThanOrEqual(numel(pk), 12, 'over-segmented');
end

function testAutoProminenceS10NoOverSegmentation(tc)
    % S10 reported 66 "peaks" under the old default (DATA_SWEEP #16) — the
    % Kα-split substrate peak counted several times.
    pk = autoPeaks(tc, ...
        'TR_S10_PTO_STO(100)_500c_150mT_20000sh_5hz_2theta omega_04162026.txt');
    verifySTOSubstrate(tc, pk);
    tc.verifyLessThanOrEqual(numel(pk), 15, 'over-segmented');
    tc.verifyTrue(all(diff(sort([pk.twoTheta])) >= 0.15), ...
        'duplicate/split picks closer than 0.15 deg');
end

function testAutoProminenceS06NoOverSegmentation(tc)
    % S06 (PTO on LAO) reported 69 "peaks" under the old default
    % (DATA_SWEEP #7). LAO substrate — only bounds asserted here.
    pk = autoPeaks(tc, ...
        'TR_S06_PTO_LAO(100)_600c_200mT_1000sh_2hz_2 theta omega_04082026.txt');
    tc.verifyGreaterThanOrEqual(numel(pk), 4);
    tc.verifyLessThanOrEqual(numel(pk), 15, 'over-segmented');
    tc.verifyTrue(all(diff(sort([pk.twoTheta])) >= 0.15), ...
        'duplicate/split picks closer than 0.15 deg');
end

function testAutoProminenceS11RealScan(tc)
    pk = autoPeaks(tc, ...
        'TR_S11_PTO_STO(100)_580c_150mT_20000sh_5hz_2theta omega_04162026.txt');
    verifySTOSubstrate(tc, pk);
    tc.verifyGreaterThanOrEqual(numel(pk), 4, 'film peaks missing');
    tc.verifyLessThanOrEqual(numel(pk), 12, 'over-segmented');
end
```

If Step 2's recorded output contradicts any bound here (e.g. a legitimate scan has 13 peaks), adjust the bound to the observed-and-judged-correct value AND say so in the task report — bounds must encode the calibrated reality, not hope. If a diff-spacing assertion conflicts with a resolved Kα₂ at STO 003 (~0.21° spacing), the 0.15° floor still passes it; do not loosen below 0.15 without evidence.

- [ ] **Step 5: Run testPeaks** — gated tests must PASS (not skip) on this machine, since `data/` is present. Run: `matlab -batch "results = runtests('tests/testPeaks.m'); disp(table(results)); assertSuccess(results)"`. Confirm in the output that the five `testAutoProminence*` tests ran (not filtered/skipped).

- [ ] **Step 6: Run the FULL suite** — same command as Task 1 Step 5. Expected: green.

- [ ] **Step 7: Commit**

```bash
git add validation/runAutoProminence.m xrdc-matlab/tests/testPeaks.m
git commit -m "test(peaks): calibrate auto prominence on Tushar data dump; gated regression tests"
```

(If Step 3 changed constants, include `xrdc-matlab/+xrdc/+peaks/findPeaks.m` in the same commit and mention the final constants in the commit body.)

---

### Task 6: Documentation

**Files:**
- Modify: `xrdc-matlab/docs/USER_GUIDE.md`
- Modify: `xrdc-matlab/docs/SCIENTIFIC_ASSUMPTIONS.md`
- Modify: `xrdc-matlab/docs/FEATURES.md`
- Modify: `validation/DATA_SWEEP.md`

- [ ] **Step 1: USER_GUIDE**

(a) In §4.6 (line ~178), change the example:

```matlab
pk   = xrdc.peaks.findPeaks(scan, 'MinProminence', 0.005*max(scan.counts));
```

to:

```matlab
pk   = xrdc.peaks.findPeaks(scan);   % auto: log-domain prominence
```

(b) After the §4.6 GUI paragraph (line ~193-195, "In the GUI, load a θ-2θ scan…"), append:

```markdown
Peak detection defaults to **auto** prominence: a peak must rise ≥ 0.3 decades
above its neighbouring troughs on the log-intensity curve and clear a 5σ Poisson
test against the local background, so a 10³-count film peak and a 10⁶-count
substrate peak are both found with no threshold tuning. The GUI's **Min prom (%)**
field reads `auto` by default; type a number to switch back to a manual threshold
(that percentage of the maximum count). The auto path also merges the Kα₁/Kα₂
substrate split (0.2° separation).
```

(c) In the §4.x θ-2θ "Core calls" block (line ~117), after the code block ending at line ~121, add the sentence:

```markdown
Omit `'MinProminence'` entirely to use the automatic log-domain criterion (see §4.6).
```

(d) Check §7 (divergence list from the Delphi original) and append a bullet:

```markdown
- **Automatic peak prominence (new, no Delphi equivalent):** when `MinProminence`
  is omitted, `findPeaks` uses a log-domain prominence criterion (≥ 0.3 decades
  above neighbouring troughs) with a 5σ Poisson noise guard, instead of a fixed
  linear threshold. The Delphi tool only had slope-threshold detectors with
  manually tuned sensitivities.
```

(Adapt wording to the list's existing bullet style after reading §7.)

- [ ] **Step 2: SCIENTIFIC_ASSUMPTIONS §3.3**

Replace the body of §3.3 (lines ~342–363: title, current-implementation block, assumption text, research prompt) with:

```markdown
### 3.3 Auto prominence: log-domain criterion (replaced the 1.5% default, 2026-06-12)

**File:** `+xrdc/+peaks/findPeaks.m` (`autoDetect`)

**Current implementation.** When `MinProminence` is omitted, a peak must
(a) rise ≥ 0.3 decades above its neighbouring troughs on `log10(max(counts,1))`
and (b) satisfy `counts − bg ≥ 5·sqrt(max(bg,1))` with `bg = movmedian(counts, ~1°)`
— i.e. a scale-invariant log-prominence test plus a Poisson significance guard.
Explicit numeric `MinProminence` keeps the fixed linear threshold.

**The assumptions.**
- 0.3 decades (~2×) reflects how peaks are judged on the standard log-intensity
  plot; calibrated on the Paik-group PTO/STO and PTO/LAO θ-2θ scans
  (validation/DATA_SWEEP.md data dump). Counting-statistics-limited features
  smaller than 2× local troughs are not reported.
- The 5σ Poisson guard assumes raw counts (not counts-per-second, not
  background-subtracted). For cps data the guard is conservative (over-rejects)
  when dwell > 1 s and permissive when dwell < 1 s.
- Counts are clamped at 1 before the log, so sub-1-count baselines produce no
  peaks; intensities below ~1 count are treated as noise by construction.
- The ~1° `movmedian` background window assumes peak FWHM ≪ 1°; correct for
  epitaxial-film HRXRD, biased for very broad (>0.5° FWHM) amorphous humps.

**Old default for reference (removed):** `max(1, 0.015·(max−min))` — missed
low-count film peaks next to 10⁶-count substrate peaks and over-segmented sharp
intense peaks (DATA_SWEEP findings, 2026).
```

Also update the summary-table row (line ~384) from
`| 3 | 1.5% prominence default | findPeaks.m:86-88 |` to
`| 3 | Log-domain auto prominence (0.3 dec / 5σ) | findPeaks.m (autoDetect) |`.

- [ ] **Step 3: FEATURES.md** — after line 34 (`- ✅ Peak detection — findPeaks …`), insert:

```markdown
- ✅ Auto prominence — log-domain (0.3 decades) + Poisson guard; default when `MinProminence` omitted
```

- [ ] **Step 4: DATA_SWEEP.md** — append to finding #3 (line ~73–76):

```markdown
   _Update 2026-06-12: fixed — `findPeaks` now defaults to a log-domain auto
   prominence with a Poisson noise guard, and the GUI θ-2θ path adds a 0.2°
   separation (Kα merge). Gated regression tests in `tests/testPeaks.m` pin
   S04/S05/S06/S10/S11 behaviour. See
   `docs/superpowers/specs/2026-06-12-auto-prominence-design.md`._
```

- [ ] **Step 5: Run the FULL suite** (docs only, but house rule) — same command as Task 1 Step 5. Expected: green.

- [ ] **Step 6: Commit**

```bash
git add xrdc-matlab/docs/USER_GUIDE.md xrdc-matlab/docs/SCIENTIFIC_ASSUMPTIONS.md xrdc-matlab/docs/FEATURES.md validation/DATA_SWEEP.md
git commit -m "docs: auto min-prominence — user guide, assumptions, features, sweep note"
```

---

## Done criteria (whole feature)

- `runtests` fully green from `xrdc-matlab/`, with the five gated `testAutoProminence*` tests RUNNING (data dump is present on this machine), not skipping.
- `xrdc.peaks.findPeaks(scan)` on S05 finds STO 001/002/003 + film peaks with zero options.
- GUI θ-2θ panel shows `auto` and plots substrate + film markers on load; typing `5` reproduces the old behaviour.
- Explicit-`MinProminence` callers (superlattice/validation scripts, XRR, RC paths) behaviourally untouched.
- Docs updated in the four files above; spec committed at `docs/superpowers/specs/2026-06-12-auto-prominence-design.md`.
