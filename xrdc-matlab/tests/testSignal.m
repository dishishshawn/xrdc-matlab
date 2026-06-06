function tests = testSignal
%TESTSIGNAL  Unit tests for the +xrdc.+signal package.
    tests = functiontests(localfunctions);
end

function testSmoothPreservesConstant(tc)
    counts = 100 * ones(50, 1);
    y = xrdc.signal.smoothCounts(counts, 5);
    tc.verifyEqual(y, counts, 'AbsTol', 1e-12);
end

function testSmoothMeanOfLinear(tc)
    % For a linear ramp, movmean preserves the middle exactly
    counts = (1:100).';
    y = xrdc.signal.smoothCounts(counts, 3);
    tc.verifyEqual(y(10:end-10), counts(10:end-10), 'AbsTol', 1e-10);
end

function testSubtractBackgroundFlatLine(tc)
    counts = 500 * ones(100, 1);
    y = xrdc.signal.subtractBackground(counts, 20);
    tc.verifyEqual(y, zeros(100, 1), 'AbsTol', 1e-10);
end

function testSubtractBackgroundClipsNegative(tc)
    counts = [10; 0; 0; 0; 10];
    y = xrdc.signal.subtractBackground(counts, 3);
    tc.verifyTrue(all(y >= 0));
end

function testDerivativesRampFirstDeriv(tc)
    % d/dx of x should be 1 everywhere (modulo edges)
    x = (0:0.01:10).';
    y = 3 * x;    % slope 3
    [slope, slope2] = xrdc.signal.derivatives(x, y, 11, 3);
    mid = 20:numel(x)-20;
    tc.verifyEqual(slope(mid),  3*ones(size(mid.')), 'AbsTol', 1e-8);
    tc.verifyEqual(slope2(mid), zeros(size(mid.')), 'AbsTol', 1e-8);
end

function testDerivativesQuadraticSecondDeriv(tc)
    % d²/dx² of x² is 2 everywhere
    x = (0:0.01:10).';
    y = x.^2;
    [~, slope2] = xrdc.signal.derivatives(x, y, 15, 3);
    mid = 50:numel(x)-50;
    tc.verifyEqual(slope2(mid), 2*ones(size(mid.')), 'AbsTol', 1e-6);
end

function testSmoothCountsCsapsPreservesLinear(tc)
    if isempty(which('csaps'))
        tc.assumeFail('csaps not available (Curve Fitting Toolbox missing).');
    end
    counts = (1:200).';
    y = xrdc.signal.smoothCounts(counts, 11, 'csaps');
    mid = 20:numel(counts)-20;
    tc.verifyEqual(y(mid), counts(mid), 'AbsTol', 1e-6);
end

function testSubtractBackgroundSplineFlatLine(tc)
    if isempty(which('csaps'))
        tc.assumeFail('csaps not available (Curve Fitting Toolbox missing).');
    end
    counts = 500 * ones(120, 1);
    y = xrdc.signal.subtractBackground(counts, 20, 'spline');
    tc.verifyEqual(y, zeros(120, 1), 'AbsTol', 1e-6);
end

function testLogEnvelopeFollowsPolynomial(tc)
    % A cubic decay should be matched by the sgolay envelope (poly order 3)
    % to within numerical noise across the interior.
    x = (0:0.01:5).';
    yLog = -2 * x + 0.1 * x.^2 - 0.01 * x.^3 + 4;
    env = xrdc.signal.logEnvelope(yLog, 41);
    mid = 100:numel(x)-100;
    tc.verifyEqual(env(mid), yLog(mid), 'AbsTol', 1e-3);
end

function testLogEnvelopeRemovesRipple(tc)
    % Slow decay + a fast sinusoidal ripple (period << envelope span):
    % the envelope must remove the decay while leaving the ripple
    % amplitude largely intact. Real XRR fringe periods are typically
    % well below the 0.4° envelope span — this mirrors that ratio.
    x = (0:0.002:3).';
    slow  = -3 * x + 5;
    ripple = 0.05 * sin(2 * pi * 25 * x);   % period ≈ 20 samples
    yLog  = slow + ripple;
    env = xrdc.signal.logEnvelope(yLog, 81);  % span ≈ 4× ripple period
    detrended = yLog - env;
    mid = 100:numel(x)-100;
    tc.verifyEqual(max(abs(detrended(mid))), 0.05, 'RelTol', 0.1);
end
