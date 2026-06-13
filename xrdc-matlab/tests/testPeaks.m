function tests = testPeaks
%TESTPEAKS  Unit tests for the +xrdc.+peaks package.
    tests = functiontests(localfunctions);
end

% ---------- helpers ----------

function scan = syntheticScan(x, centers, fwhms, amps, shape, bgLevel, noiseSd)
    if nargin < 6, bgLevel = 100;  end
    if nargin < 7, noiseSd = 0;    end
    x = x(:);
    y = bgLevel * ones(size(x));
    for k = 1:numel(centers)
        y = y + evalShape(x, centers(k), fwhms(k), amps(k), shape);
    end
    if noiseSd > 0
        rng(1);
        y = y + noiseSd * randn(size(y));
    end
    scan = xrdc.io.emptyScan();
    scan.twoTheta = x;
    scan.counts   = y;
    scan.sourceFormat = "synthetic";
end

function y = evalShape(x, x0, fwhm, amp, shape)
    switch shape
        case "lorentz"
            g = fwhm / 2;
            y = amp * g^2 ./ ((x - x0).^2 + g^2);
        case "gauss"
            sigma = fwhm / (2*sqrt(2*log(2)));
            y = amp * exp(-(x - x0).^2 / (2*sigma^2));
    end
end

% ---------- findPeaks (modern) ----------

function testFindPeaksSingleLorentz(tc)
    x = (20:0.02:40).';
    scan = syntheticScan(x, 30, 0.3, 10000, "lorentz", 200);
    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyLength(pk, 1);
    tc.verifyEqual(pk(1).twoTheta, 30, 'AbsTol', 0.05);
    tc.verifyEqual(pk(1).fwhm,     0.3, 'RelTol', 0.2);
end

function testFindPeaksMultiplePeaks(tc)
    x = (20:0.02:80).';
    scan = syntheticScan(x, [28, 47, 69], [0.3, 0.25, 0.35], ...
                              [8000, 6000, 10000], "gauss", 150);
    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyGreaterThanOrEqual(numel(pk), 3);
    % Check that the 3 strongest are at the expected centres
    [~, order] = sort([pk.counts], 'descend');
    top = sort([pk(order(1:3)).twoTheta]);
    tc.verifyEqual(top, [28, 47, 69], 'AbsTol', 0.05);
end

function testFindPeaksHeightThreshold(tc)
    x = (20:0.02:40).';
    scan = syntheticScan(x, [25, 35], [0.3, 0.3], [100, 5000], ...
                              "lorentz", 50);
    pk = xrdc.peaks.findPeaks(scan, 'MinHeight', 1000);
    tc.verifyLength(pk, 1);
    tc.verifyEqual(pk(1).twoTheta, 35, 'AbsTol', 0.05);
end

function testFindPeaksRangeCrop(tc)
    x = (20:0.02:80).';
    scan = syntheticScan(x, [28, 47, 69], [0.3, 0.25, 0.35], ...
                              [8000, 6000, 10000], "gauss", 150);
    pk = xrdc.peaks.findPeaks(scan, 'TwoThetaRange', [40, 60]);
    tc.verifyLength(pk, 1);
    tc.verifyEqual(pk(1).twoTheta, 47, 'AbsTol', 0.05);
end

function testFindPeaksEmptyScan(tc)
    scan = xrdc.io.emptyScan();
    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyEmpty(pk);
end

function testFindPeaksBadScan(tc)
    scan = struct('counts', [1;2;3]);   % no twoTheta field
    tc.verifyError(@() xrdc.peaks.findPeaks(scan), ...
        'xrdc:peaks:badScan');
end

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

function testFindPeaksAutoRejectsQuantisedCpsBlips(tc)
    % Regression for the cps-units / mostly-zero-baseline bug.
    %
    % Real Rigaku θ-2θ exports are counts-per-second: one detected photon
    % is ~5 cps, not 1 count, and ~70% of the baseline is exactly 0. The
    % old guard assumed σ = sqrt(counts) on the raw cps and floored σ at 1
    % count, so a single-photon (5 cps) baseline blip cleared a flat 5σ
    % threshold and was reported as a peak — ~647 of them on the S05 scan.
    %
    % Here we emulate that: counts are mostly 0 with isolated single-quantum
    % (q = 5) blips, PLUS two genuine peaks — a ~1e6 substrate-like and a
    % ~1e3 film-like Gaussian. Quantise to multiples of q so the data looks
    % like real cps (the quantum estimator must recover q ≈ 5). Auto mode
    % must return exactly the two real peaks and none of the blips.
    q = 5;
    x = (20:0.01:50).';
    n = numel(x);
    y = 1e6 * exp(-(x - 22.7).^2 / (2 * 0.02^2)) ...   % substrate-like
      + 1e3 * exp(-(x - 46.5).^2 / (2 * 0.10^2));       % film-like
    rng(12345);
    nBlip = 150;
    blipIdx = randperm(n, nBlip);
    y(blipIdx) = y(blipIdx) + q;                          % single-quantum blips
    y = q * round(y / q);                                 % quantise like cps
    scan = xrdc.io.emptyScan();
    scan.twoTheta = x;
    scan.counts   = y;
    scan.sourceFormat = "synthetic";

    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyLength(pk, 2);
    tc.verifyEqual(sort([pk.twoTheta]), [22.7, 46.5], 'AbsTol', 0.05);
end

function testFindPeaksAutoMarginalLogProminencePeak(tc)
    % Constant-pin PROM_DECADES: a peak only ~0.4 decades above its
    % neighbouring troughs (just clear of the 0.3-decade shape test) must
    % still be detected. If PROM_DECADES ever drifts up past ~0.4 this fails.
    % Background 1000, peak amplitude 1512 -> top ~2512 ≈ 0.4 decades above
    % the 1000 troughs (log10(2512/1000) ≈ 0.40). Strong Poisson
    % significance keeps the noise guard from being the limiting factor.
    x = (20:0.01:40).';
    scan = syntheticScan(x, 30, 0.4, 1512, "gauss", 1000);
    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyGreaterThanOrEqual(numel(pk), 1);
    [~, ord] = sort([pk.counts], 'descend');
    tc.verifyEqual(pk(ord(1)).twoTheta, 30, 'AbsTol', 0.05);
end

function testFindPeaksAutoWeakPeakAboveModestBackground(tc)
    % Constant-pin NOISE_SIGMAS: a weak-but-real peak must survive the
    % Poisson significance guard. Background 50, narrow (FWHM 0.2°) peak of
    % amplitude 75 -> top 125. The log shape test is comfortably cleared
    % (log10(125/50) ≈ 0.40 decades, well above 0.3), so the binding
    % constraint is the noise guard alone. The photon-excess significance
    % here is z = (125 - 50)/sqrt(125 + 50 + 1) ≈ 5.6 (q = 1 on this
    % continuous synthetic data), just above the 5σ guard. If NOISE_SIGMAS
    % ever drifts up past ~5.5 this genuine peak is rejected and the test
    % fails — pinning the guard from becoming too strict (well under the
    % review's ≤ ~8 bound).
    x = (20:0.01:40).';
    scan = syntheticScan(x, 30, 0.2, 75, "gauss", 50);
    pk = xrdc.peaks.findPeaks(scan);
    tc.verifyGreaterThanOrEqual(numel(pk), 1);
    [~, ord] = sort([pk.counts], 'descend');
    tc.verifyEqual(pk(ord(1)).twoTheta, 30, 'AbsTol', 0.05);
end

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

% ---------- findPeaksLegacy (slope / slope2) ----------

function testFindPeaksLegacySlopeSingle(tc)
    x = (20:0.02:40).';
    scan = syntheticScan(x, 30, 0.3, 10000, "lorentz", 200);
    pk = xrdc.peaks.findPeaksLegacy(scan, 'SlopeThreshold', 500);
    tc.verifyGreaterThanOrEqual(numel(pk), 1);
    tt = [pk.twoTheta];
    tc.verifyTrue(any(abs(tt - 30) < 0.3), ...
        'Legacy slope detector missed the only peak at 30°.');
end

function testFindPeaksLegacySlope2Single(tc)
    x = (20:0.02:40).';
    scan = syntheticScan(x, 30, 0.3, 10000, "lorentz", 200);
    pk = xrdc.peaks.findPeaksLegacy(scan, ...
        'Algorithm', 'slope2', 'SlopeThreshold', 100);
    tc.verifyGreaterThanOrEqual(numel(pk), 1);
    tt = [pk.twoTheta];
    tc.verifyTrue(any(abs(tt - 30) < 0.3));
end

function testFindPeaksLegacyMinSeparation(tc)
    % Two peaks 0.02° apart — within any reasonable merge tolerance
    x = (20:0.005:40).';
    scan = syntheticScan(x, [30, 30.02], [0.1, 0.1], [5000, 5000], ...
                              "lorentz", 100);
    pk = xrdc.peaks.findPeaksLegacy(scan, ...
        'SlopeThreshold', 500, 'MinSeparation', 0.1);
    % Should collapse to one peak
    tc.verifyLessThanOrEqual(numel(pk), 2);
end

% ---------- adjustPeaks ----------

function testAdjustPeaksRecoversCenter(tc)
    % Deliberately mis-identify the peak centre by one step and verify
    % adjustPeaks snaps it back to 0.01° of truth.
    x = (20:0.02:40).';
    scan = syntheticScan(x, 30, 0.4, 8000, "gauss", 100);
    approx = struct('twoTheta', 30.08, 'counts', interp1(x, scan.counts, 30.08), ...
                    'prominence', NaN, 'fwhm', NaN, 'leftHalf', NaN, ...
                    'rightHalf', NaN, 'index', NaN);
    refined = xrdc.peaks.adjustPeaks(scan, approx);
    tc.verifyEqual(refined.twoTheta, 30, 'AbsTol', 0.02);
    tc.verifyEqual(refined.fwhm,    0.4, 'RelTol', 0.15);
    tc.verifyGreaterThan(refined.fwhm, 0);
end

function testAdjustPeaksMergesDuplicates(tc)
    x = (20:0.02:40).';
    scan = syntheticScan(x, 30, 0.3, 5000, "lorentz", 100);
    p = emptyPkStruct();
    p1 = p; p1.twoTheta = 30.00; p1.counts = 5100;
    p2 = p; p2.twoTheta = 30.01; p2.counts = 5050;
    refined = xrdc.peaks.adjustPeaks(scan, [p1; p2], ...
        'MergeTolerance', 0.03);
    tc.verifyLength(refined, 1);
end

function testAdjustPeaksEmpty(tc)
    x = (20:0.02:40).';
    scan = syntheticScan(x, 30, 0.3, 5000, "lorentz", 100);
    refined = xrdc.peaks.adjustPeaks(scan, struct([]));
    tc.verifyEmpty(refined);
end

% ---------- fitPeak ----------

function testFitPeakLorentzRecovery(tc)
    if isempty(which('lsqcurvefit'))
        tc.assumeFail('lsqcurvefit not available (Optimization Toolbox).');
    end
    x = (28:0.01:32).';
    scan = syntheticScan(x, 30, 0.35, 20000, "lorentz", 200);
    r = xrdc.peaks.fitPeak(scan, [28, 32], 'Shape', 'lorentz');
    tc.verifyEqual(r.twoTheta,  30,    'AbsTol', 1e-3);
    tc.verifyEqual(r.fwhm,      0.35,  'RelTol', 1e-2);
    tc.verifyEqual(r.amplitude, 20000, 'RelTol', 1e-2);
    tc.verifyGreaterThan(r.rSquared, 0.999);
end

function testFitPeakGaussRecovery(tc)
    if isempty(which('lsqcurvefit'))
        tc.assumeFail('lsqcurvefit not available (Optimization Toolbox).');
    end
    x = (28:0.01:32).';
    scan = syntheticScan(x, 30, 0.25, 15000, "gauss", 150);
    r = xrdc.peaks.fitPeak(scan, [28, 32], 'Shape', 'gauss');
    tc.verifyEqual(r.twoTheta,  30,    'AbsTol', 1e-3);
    tc.verifyEqual(r.fwhm,      0.25,  'RelTol', 1e-2);
    tc.verifyEqual(r.amplitude, 15000, 'RelTol', 1e-2);
end

function testFitPeakPseudoVoigtRecovery(tc)
    if isempty(which('lsqcurvefit'))
        tc.assumeFail('lsqcurvefit not available (Optimization Toolbox).');
    end
    % Build a synthetic pseudo-Voigt (50/50 mix) and recover η
    x = (28:0.01:32).';
    fwhm = 0.3;
    yL = evalShape(x, 30, fwhm, 1, "lorentz");
    yG = evalShape(x, 30, fwhm, 1, "gauss");
    y  = 100 + 10000 * (0.5 * yL + 0.5 * yG);
    scan = xrdc.io.emptyScan();
    scan.twoTheta = x; scan.counts = y;
    r = xrdc.peaks.fitPeak(scan, [28, 32], 'Shape', 'pseudoVoigt');
    tc.verifyEqual(r.twoTheta, 30,   'AbsTol', 1e-3);
    tc.verifyEqual(r.fwhm,     fwhm, 'RelTol', 1e-2);
    tc.verifyEqual(r.eta,      0.5,  'AbsTol', 0.1);
end

function testFitPeakReportsFiniteSE(tc)
    if isempty(which('lsqcurvefit'))
        tc.assumeFail('lsqcurvefit not available (Optimization Toolbox).');
    end
    % With Gaussian noise the SEs should be finite and positive.
    rng(3);
    x = (28:0.01:32).';
    scan = syntheticScan(x, 30, 0.3, 8000, "lorentz", 200, 20);
    r = xrdc.peaks.fitPeak(scan, [28, 32], 'Shape', 'lorentz');
    tc.verifyTrue(isfinite(r.paramSE.center) && r.paramSE.center > 0);
    tc.verifyTrue(isfinite(r.paramSE.fwhm)   && r.paramSE.fwhm   > 0);
    tc.verifyTrue(isfinite(r.paramSE.amplitude) && r.paramSE.amplitude > 0);
end

function testFitPeakMultiStartRecovery(tc)
    if isempty(which('lsqcurvefit'))
        tc.assumeFail('lsqcurvefit not available (Optimization Toolbox).');
    end
    if isempty(which('MultiStart'))
        tc.assumeFail('MultiStart not available (Global Optimization Toolbox).');
    end
    rng(7);
    x = (28:0.01:32).';
    scan = syntheticScan(x, 30, 0.4, 12000, "lorentz", 200, 30);
    r = xrdc.peaks.fitPeak(scan, [28, 32], ...
        'Shape', 'lorentz', 'Method', 'multistart', 'StartPoints', 8);
    tc.verifyEqual(r.twoTheta, 30,   'AbsTol', 0.01);
    tc.verifyEqual(r.fwhm,     0.4,  'RelTol', 0.1);
    tc.verifyGreaterThan(r.rSquared, 0.95);
end

function testFitPeakBadWindow(tc)
    scan = xrdc.io.emptyScan();
    scan.twoTheta = (20:0.02:40).';
    scan.counts   = ones(size(scan.twoTheta));
    tc.verifyError(@() xrdc.peaks.fitPeak(scan, [30, 20]), ...
        'xrdc:peaks:badWindow');
end

function testFitPeakTooFewPoints(tc)
    scan = xrdc.io.emptyScan();
    scan.twoTheta = [29.98; 30.00; 30.02];   % only 3 points
    scan.counts   = [100;   200;   100];
    tc.verifyError( ...
        @() xrdc.peaks.fitPeak(scan, [29.97, 30.03]), ...
        'xrdc:peaks:tooFewPoints');
end

% ---------- findPhiPeaks ----------

function testFindPhiPeaksRejectsNoiseAndCountsFold(tc)
    % Weak 4-fold pole scan: 4 poles ~20 counts on a 1-count background,
    % plus single-bin noise spikes of 3 counts. The noise-floor threshold
    % must keep the 4 poles and reject the spikes (the legacy "10% of max"
    % rule would have counted the spikes).
    phi = (-180:0.2:179.8).';
    y   = ones(size(phi));
    for c = [-135 -45 45 135]
        y = y + 20 * exp(-((phi - c).^2) / (2 * 0.3^2));
    end
    y([100 500 900 1300]) = 3;          % noise spikes, below bg + 6σ

    s = xrdc.io.emptyScan();
    s.twoTheta = phi; s.counts = y; s.scanType = "phi";

    [pk, info] = xrdc.peaks.findPhiPeaks(s);
    tc.verifyEqual(numel(pk), 4);
    tc.verifyEqual(info.nUnique, 4);
    tc.verifyEqual(info.fold, 4);
    tc.verifyLessThan(max(abs(info.spacings - 90)), 1);
end

function testFindPhiPeaksDedupesWrapAcross360(tc)
    % An over-rotated scan (>360°) captures the same pole at both ends; it
    % should be counted once, still 4-fold.
    phi = (-200:0.2:200).';             % 400° span
    y   = ones(size(phi));
    for c = [-180 -90 0 90 180]         % +180 and -180 are the same pole
        y = y + 30 * exp(-((phi - c).^2) / (2 * 0.3^2));
    end
    s = xrdc.io.emptyScan();
    s.twoTheta = phi; s.counts = y; s.scanType = "phi";

    [~, info] = xrdc.peaks.findPhiPeaks(s);
    tc.verifyEqual(info.nUnique, 4);
    tc.verifyEqual(info.fold, 4);
end

function testFindPhiPeaksFlagsBoundaryWrap(tc)
    % Mimic the real PtO2 geometry: 4 poles 90° apart, scan spanning ~360°
    % so one pole's next-period instance has its apex just past the end,
    % leaving a cut-off edge at the boundary. That edge must be flagged as a
    % wrap repeat (same direction as a counted pole), NOT a 5th pole.
    phi = (-197.9:0.2:161.9).';
    y   = ones(size(phi));
    for c = [-197.4 -107.4 -17.4 72.6 162.6]   % 162.6 ≡ -197.4 (360° away)
        y = y + 30 * exp(-((phi - c).^2) / (2 * 0.5^2));
    end
    s = xrdc.io.emptyScan();
    s.twoTheta = phi; s.counts = y; s.scanType = "phi";

    [pk, info] = xrdc.peaks.findPhiPeaks(s);
    tc.verifyEqual(numel(pk), 4);                       % 4 poles, not 5
    tc.verifyEqual(info.fold, 4);
    tc.verifyGreaterThanOrEqual(numel(info.wrapRepeats), 1);
    % the wrap sits near the scan end, not on a counted pole
    tc.verifyGreaterThan(info.wrapRepeats(1).twoTheta, 150);
end

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
%Helper (not a test): the STO 001/002/003 triplet present, once each, at nominal.
    tt = [pk.twoTheta];
    for ref = [22.75, 46.47, 72.57]
        tc.verifyEqual(nnz(abs(tt - ref) < 0.05), 1, ...
            sprintf('expected exactly one peak within 0.05 deg of %.2f', ref));
    end
end

function testAutoProminenceS05RealScan(tc)
    % S05: previously needed a manually lowered threshold. Conservative
    % auto returns the STO substrate triplet (film below noise).
    pk = autoPeaks(tc, ...
        'TR_S05_PTO_STO(100)_600c_200mT_1000sh_2hz_2theta omega_04092026.txt');
    verifySTOSubstrate(tc, pk);
    tc.verifyGreaterThanOrEqual(numel(pk), 3);
    tc.verifyLessThanOrEqual(numel(pk), 12, 'over-segmented');
    tc.verifyTrue(all(diff(sort([pk.twoTheta])) >= 0.15), ...
        'duplicate/split picks closer than 0.15 deg');
end

function testAutoProminenceS04RealScan(tc)
    pk = autoPeaks(tc, ...
        'TR_S04_PTO_STO(100)_750c_200mT_1000sh_3hz_2theta omega_04072026.txt');
    verifySTOSubstrate(tc, pk);
    tc.verifyGreaterThanOrEqual(numel(pk), 3);
    tc.verifyLessThanOrEqual(numel(pk), 12, 'over-segmented');
end

function testAutoProminenceS10NoOverSegmentation(tc)
    % S10 reported 66 "peaks" under the old default (DATA_SWEEP #16) — the
    % Kα-split substrate peak counted several times. Its substrate peaks
    % are sample-shifted ~0.1 deg off nominal, so assert the
    % over-segmentation bounds, not the 0.05-deg triplet.
    pk = autoPeaks(tc, ...
        'TR_S10_PTO_STO(100)_500c_150mT_20000sh_5hz_2theta omega_04162026.txt');
    tc.verifyGreaterThanOrEqual(numel(pk), 3);
    tc.verifyLessThanOrEqual(numel(pk), 15, 'over-segmented');
    tc.verifyTrue(all(diff(sort([pk.twoTheta])) >= 0.15), ...
        'duplicate/split picks closer than 0.15 deg');
    % substrate peaks present near (sample-shifted) nominal positions
    tt = sort([pk.twoTheta]);
    tc.verifyTrue(any(abs(tt - 46.5) < 0.2), 'STO 002 region peak missing');
end

function testAutoProminenceS06NoOverSegmentation(tc)
    % S06 (PTO on LAO) reported 69 "peaks" under the old default
    % (DATA_SWEEP #7). LAO substrate — only bounds asserted.
    pk = autoPeaks(tc, ...
        'TR_S06_PTO_LAO(100)_600c_200mT_1000sh_2hz_2 theta omega_04082026.txt');
    tc.verifyGreaterThanOrEqual(numel(pk), 2);
    tc.verifyLessThanOrEqual(numel(pk), 15, 'over-segmented');
    tc.verifyTrue(all(diff(sort([pk.twoTheta])) >= 0.15), ...
        'duplicate/split picks closer than 0.15 deg');
end

function testAutoProminenceS11FindsFilm(tc)
    % S11 has genuine strong film peaks (~10x local background) that the
    % conservative auto path DOES detect: STO triplet + PTO film.
    pk = autoPeaks(tc, ...
        'TR_S11_PTO_STO(100)_580c_150mT_20000sh_5hz_2theta omega_04162026.txt');
    verifySTOSubstrate(tc, pk);
    tc.verifyGreaterThanOrEqual(numel(pk), 5, 'film peaks missing');
    tc.verifyLessThanOrEqual(numel(pk), 12, 'over-segmented');
    tt = sort([pk.twoTheta]);
    tc.verifyTrue(any(abs(tt - 21.6) < 0.3), 'PTO film ~21.6 deg missing');
end

% ---------- small helpers ----------

function s = emptyPkStruct()
    s = struct( ...
        'twoTheta',   NaN, ...
        'counts',     NaN, ...
        'prominence', NaN, ...
        'fwhm',       NaN, ...
        'leftHalf',   NaN, ...
        'rightHalf',  NaN, ...
        'index',      NaN);
end
