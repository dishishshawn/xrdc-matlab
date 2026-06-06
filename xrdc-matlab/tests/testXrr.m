function tests = testXrr
%TESTXRR  Unit tests for the +xrdc.+xrr package.
    tests = functiontests(localfunctions);
end

% =========================================================================
% Helpers
% =========================================================================

function scan = syntheticXrr(tNm, lambda, twoThetaC, twoThetaEnd)
%SYNTHETICXRR  Synthetic Rigaku-style XRR with direct beam, TER plateau
%   peaking at twoThetaC, post-edge decay, and Kiessig fringes for a
%   film of thickness tNm (nm).
    if nargin < 3, twoThetaC   = 0.80; end
    if nargin < 4, twoThetaEnd = 5.0;  end
    step = 0.004;
    tt   = (0:step:twoThetaEnd).';
    n    = numel(tt);

    % Direct-beam Gaussian spike at 2θ = 0 (width ~0.02°)
    directBeam = 3e6 * exp(-(tt / 0.02).^2);

    % Smooth envelope peaking at twoThetaC, then decaying as (θ_c/θ)^4
    rise = exp(-((tt - twoThetaC) / 0.10).^2);
    decay = ones(n, 1);
    post = tt > twoThetaC;
    decay(post) = (twoThetaC ./ tt(post)).^4;
    envelope = 5e5 * max(rise, decay);

    % Kiessig fringes only kick in past the critical edge
    lambdaNm = lambda / 10;
    sinTh = sin(tt * pi / 180 / 2);
    phase = (sinTh - sin(twoThetaC * pi / 180 / 2)) * 2 * tNm / lambdaNm;
    fringes = zeros(n, 1);
    fringes(post) = 0.4 * cos(2 * pi * phase(post));

    counts = directBeam + envelope .* (1 + fringes);
    counts = max(counts, 1);

    scan = xrdc.io.emptyScan();
    scan.twoTheta     = tt;
    scan.counts       = counts;
    scan.lambda       = lambda;
    scan.scanType     = "twoThetaOmega";
    scan.identifier   = sprintf("synthetic XRR t=%g nm", tNm);
    scan.sourceFormat = "synthetic";
end

% =========================================================================
% findCriticalEdge
% =========================================================================

function testCriticalEdgeOnSynthetic(tc)
    scan = syntheticXrr(30, 1.5406, 0.80);
    info = xrdc.xrr.findCriticalEdge(scan.twoTheta, scan.counts);

    % Direct-beam dip should be in the very first 0.2°
    tc.verifyLessThan(info.twoThetaDbDip, 0.2);
    % TER plateau peak should be at the synthetic θ_c (0.80°).
    tc.verifyEqual(info.twoThetaC, 0.80, 'AbsTol', 0.10);
end

function testCriticalEdgeRejectsTinyScan(tc)
    twoTheta = (0:0.1:0.3).';
    counts   = ones(size(twoTheta));
    tc.verifyError( ...
        @() xrdc.xrr.findCriticalEdge(twoTheta, counts), ...
        'xrdc:xrr:tooFewPoints');
end

function testCriticalEdgeRejectsMismatchedInputs(tc)
    tc.verifyError( ...
        @() xrdc.xrr.findCriticalEdge((1:10).', (1:11).'), ...
        'xrdc:xrr:sizeMismatch');
end

% =========================================================================
% detrendLog
% =========================================================================

function testDetrendLogRemovesPolynomial(tc)
    x = linspace(0.5, 4, 500).';
    y = 5 - 2 * x + 0.3 * x.^2 - 0.02 * x.^3;     % a 4th-order polynomial
    yDet = xrdc.xrr.detrendLog(x, y, 'Order', 4);
    tc.verifyEqual(yDet, zeros(size(y)), 'AbsTol', 1e-6);
end

function testDetrendLogPreservesRipple(tc)
    x  = linspace(0.5, 4, 1000).';
    slow   = 5 - 2 * x;
    ripple = 0.05 * sin(2 * pi * 8 * x);          % 8 cycles/° = 70 nm fringe
    yDet   = xrdc.xrr.detrendLog(x, slow + ripple, 'Order', 4);
    mid    = 50:numel(x)-50;
    tc.verifyEqual(max(abs(yDet(mid))), 0.05, 'RelTol', 0.10);
end

% =========================================================================
% dominantPeriod
% =========================================================================

function testDominantPeriodRecoversSinusoid(tc)
    step = 0.004;
    x    = (0:step:4).';
    f    = 5;                                      % cycles per degree
    y    = sin(2 * pi * f * x);
    info = xrdc.xrr.dominantPeriod(x, y);
    tc.verifyEqual(info.freqPerDeg, f, 'RelTol', 0.05);
    tc.verifyEqual(info.periodDeg,  1/f, 'RelTol', 0.05);
end

function testDominantPeriodRejectsTinyVector(tc)
    tc.verifyError( ...
        @() xrdc.xrr.dominantPeriod((1:5).', (1:5).'), ...
        'xrdc:xrr:tooFewPoints');
end

% =========================================================================
% analyzeFringes (end-to-end on synthetic data)
% =========================================================================

function testAnalyzeFringesSyntheticThickness(tc)
    % Build a 30 nm synthetic film and confirm the pipeline recovers it
    % within a few percent on both the FFT and the linear-fit thickness.
    scan = syntheticXrr(30, 1.5406);
    res  = xrdc.xrr.analyzeFringes(scan);

    tc.verifyEqual(res.thicknessFftNm, 30, 'RelTol', 0.10);
    tc.verifyGreaterThanOrEqual(res.nFringesDetected, 5);
    tc.verifyEqual(res.thicknessFitNm, 30, 'RelTol', 0.10);
end

function testAnalyzeFringesManualLowerBound(tc)
    scan = syntheticXrr(30, 1.5406);
    res  = xrdc.xrr.analyzeFringes(scan, 'LowerBound', 1.0);
    tc.verifyEqual(res.lowerBound, 1.0, 'AbsTol', 1e-9);
end

function testThicknessQuadraticRecoversSyntheticThickness(tc)
    scan = syntheticXrr(35, 1.5406, 0.80);
    res  = xrdc.xrr.analyzeFringes(scan);
    tc.verifyEqual(res.thicknessQuadNm, 35, 'RelTol', 0.05);
    % quadratic should agree with FFT within a couple of percent
    tc.verifyEqual(res.thicknessQuadNm, res.thicknessFftNm, 'RelTol', 0.05);
end

function testThicknessQuadraticRequires3Fringes(tc)
    tc.verifyError( ...
        @() xrdc.xrr.thicknessQuadratic([0.7; 1.0], 1.5406), ...
        'xrdc:xrr:tooFewFringes');
end

function testAnalyzeFringesRealSampleS25(tc)
    % Regression on the user's S25 sample. Pre-fix the pipeline returned
    % 105.5 nm. The XRR data themselves — under three independent methods
    % (linear sinθ, quadratic Kiessig, FFT) — converge on ~40 nm with
    % residuals 1e-5, so the regression target is 40 nm, not the nominal
    % 36 nm from PLD calibration. PLD growth-rate calibration commonly
    % carries ±10% error; the XRR fringe analysis is the more accurate
    % measurement.
    candidates = { ...
        ['C:\Users\TNTMi\Desktop\XDRC MATLAB Project\tushar_data\' ...
         'TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_XRR_05062026.txt'], ...
        ['C:\Users\TNTMi\Desktop\XDRC MATLAB Project\New folder\' ...
         'TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_XRR_05062026.txt']};
    fpath = '';
    for k = 1:numel(candidates)
        if isfile(candidates{k}), fpath = candidates{k}; break; end
    end
    tc.assumeFalse(isempty(fpath), ...
        'S25 reference data not present on this machine.');

    scan = xrdc.io.readScan(fpath);
    res  = xrdc.xrr.analyzeFringes(scan);

    tc.verifyEqual(res.thicknessQuadNm, 40, 'AbsTol', 1.5);
    tc.verifyEqual(res.thicknessFftNm,  40, 'AbsTol', 1.5);
    tc.verifyGreaterThan(res.fringeSnr, 3);
end
