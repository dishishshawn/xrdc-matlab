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
        y = y + evalShapeTest(x, centers(k), fwhms(k), amps(k), shape);
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

function y = evalShapeTest(x, x0, fwhm, amp, shape)
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
    yL = evalShapeTest(x, 30, fwhm, 1, "lorentz");
    yG = evalShapeTest(x, 30, fwhm, 1, "gauss");
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
