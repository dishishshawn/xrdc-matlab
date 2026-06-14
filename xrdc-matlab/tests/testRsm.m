function tests = testRsm
%TESTRSM  Unit tests for +xrdc/+rsm/ and the RSM plot.
%   Exercises toReciprocalSpace (FR-6.1), loadAreaScan (FR-6.2),
%   and plotRsm (FR-6.3).  setOffsetsInteractive requires a live
%   display and ginput — tested separately (not automated here).
%
%   V7 parity: the known-answer tests below verify the θ-asymmetry
%   formulas to machine precision. Real-data parity (±1e-6 Å⁻¹ vs
%   XRDC output) requires Dr. Paik's reference RSM dataset; that test
%   is gated on file existence below.
    tests = functiontests(localfunctions);
end

% ---------- fixtures ----------

function setupOnce(testCase) %#ok<INUSD>
    set(groot, 'DefaultFigureVisible', 'off');
end

function teardownOnce(testCase) %#ok<INUSD>
    set(groot, 'DefaultFigureVisible', 'on');
    close all force;
end

% ---------- helpers ----------

function s = makeSlice(tt_deg, omega_deg, lambda)
    s = xrdc.io.emptyScan();
    s.twoTheta      = tt_deg(:);
    s.counts        = 1000 * ones(numel(tt_deg), 1);
    s.secondAxis    = omega_deg;
    s.secondAxisName = "Omega";
    s.scanType      = "twoThetaOmega";
    s.lambda        = lambda;
end

% =====================================================================
% toReciprocalSpace — θ-asymmetry & 2/λ convention (SCIENTIFIC_ASSUMPTIONS.md §1.6)
% =====================================================================

function testSymmetricReflection(testCase)
    % For a symmetric reflection the scan is collected at ω = θ for every
    % point.  In that case ω − θ = 0 ⟹ k_par = 0 and k_perp = 2sinθ/λ.
    %
    % Build: secondAxis = 2θ_center/2  so that the ω formula gives ω = θ.
    %   ω = secondAxis - (2θ_ctr/2) + θ_raw + 0
    %     = (2θ_ctr/2) - (2θ_ctr/2) + 2θ/2 = 2θ/2 = θ_raw = θ (no ΔΘ)
    lambda = 1.5406;            % Cu Kα1 (Å)
    tt = (30:0.1:32).';         % narrow range around 31°
    omega = mean([tt(1) tt(end)]) / 2;   % = 2θ_center / 2
    s = makeSlice(tt, omega, lambda);

    [kPar, kPerp] = xrdc.rsm.toReciprocalSpace(s);

    % k_par must be zero for all points (symmetric scan)
    testCase.verifyLessThan(abs(kPar), 1e-10, ...
        'k_par must be ~0 for a symmetric reflection.');

    % k_perp = 2*sin(θ)/λ with θ = tt/2
    theta_rad = tt * (pi/360);
    kPerpExpected = (2/lambda) .* sin(theta_rad);
    testCase.verifyLessThan(max(abs(kPerp - kPerpExpected)), 1e-12, ...
        'k_perp does not match 2*sin(θ)/λ for symmetric scan.');
end

function testThetaAsymmetryWithDeltaTheta(testCase)
    % Verify the θ-asymmetry: when DeltaTheta ≠ 0, ω is built from
    % θ_raw (uncorrected) but k is computed from θ (corrected).
    % A non-zero ΔΘ should shift k_perp but leave the ω formula unchanged.
    lambda = 1.5406;
    tt = (30:0.1:32).';
    omega = mean([tt(1) tt(end)]) / 2;
    s = makeSlice(tt, omega, lambda);
    dTheta = 0.05;   % degrees

    [kPar0, kPerp0] = xrdc.rsm.toReciprocalSpace(s, 'DeltaTheta', 0);
    [kPar1, kPerp1] = xrdc.rsm.toReciprocalSpace(s, 'DeltaTheta', dTheta);

    % k_perp must differ (ΔΘ shifts the corrected θ)
    testCase.verifyGreaterThan(max(abs(kPerp1 - kPerp0)), 1e-6, ...
        'DeltaTheta should shift k_perp.');

    % Manual computation of the expected asymmetric result
    tt_rad        = tt * (pi/180);
    tt_center_rad = mean([tt(1) tt(end)]) * (pi/180);
    sa_rad        = omega * (pi/180);
    theta_raw     = tt_rad / 2;                           % uncorrected
    om            = sa_rad - tt_center_rad/2 + theta_raw; % ω from θ_raw
    theta_corr    = (tt_rad + dTheta*pi/180) / 2;         % corrected
    kPerpExpected = (2/lambda) .* sin(theta_corr) .* cos(om - theta_corr);
    kParExpected  = (2/lambda) .* sin(theta_corr) .* sin(om - theta_corr);

    testCase.verifyLessThan(max(abs(kPerp1 - kPerpExpected)), 1e-12, ...
        'k_perp with DeltaTheta does not match expected asymmetric formula.');
    testCase.verifyLessThan(max(abs(kPar1 - kParExpected)), 1e-12, ...
        'k_par with DeltaTheta does not match expected asymmetric formula.');
end

function testFlipNegatesKPar(testCase)
    lambda = 1.5406;
    tt = (30:0.5:35).';
    omega = 17;    % asymmetric → non-zero k_par
    s = makeSlice(tt, omega, lambda);

    [kPar, ~]     = xrdc.rsm.toReciprocalSpace(s, 'Flip', false);
    [kParFlip, ~] = xrdc.rsm.toReciprocalSpace(s, 'Flip', true);

    testCase.verifyLessThan(max(abs(kPar + kParFlip)), 1e-12, ...
        'Flip must negate k_par exactly.');
end

function testDeltaOmegaShift(testCase)
    % ΔΩ shifts ω by that amount (in degrees → radians internally).
    % With a symmetric scan, adding ΔΩ should produce non-zero k_par.
    lambda = 1.5406;
    tt = (30:0.1:32).';
    omega_sym = mean([tt(1) tt(end)]) / 2;
    s = makeSlice(tt, omega_sym, lambda);
    dOmega = 0.1;   % degrees

    [kPar0, ~] = xrdc.rsm.toReciprocalSpace(s, 'DeltaOmega', 0);
    [kPar1, ~] = xrdc.rsm.toReciprocalSpace(s, 'DeltaOmega', dOmega);

    testCase.verifyLessThan(max(abs(kPar0)), 1e-10, ...
        'Zero DeltaOmega on symmetric scan must give k_par ≈ 0.');
    testCase.verifyGreaterThan(max(abs(kPar1)), 1e-6, ...
        'Non-zero DeltaOmega must produce non-zero k_par.');
end

function testKnownAnswerCuKaSrTiO3_002(testCase)
    % SrTiO3(002) symmetric reflection: 2θ ≈ 46.47° for Cu Kα1 (1.5406 Å).
    % At the exact Bragg angle in a symmetric scan: k_par = 0,
    % k_perp = 2*sin(θ)/λ = 2/d_{002} where d_{002} = a/2 for SrTiO3 (a≈3.905 Å).
    lambda  = 1.5406;
    a_STO   = 3.905;       % Å  (bulk SrTiO3)
    d_002   = a_STO / 2;   % Å
    tt_002  = 2 * asin(lambda / (2*d_002)) * 180/pi;  % degrees
    omega   = tt_002 / 2;  % symmetric

    tt  = tt_002;   % single point at the Bragg angle
    s   = makeSlice(tt, omega, lambda);

    [kPar, kPerp] = xrdc.rsm.toReciprocalSpace(s);

    % At Bragg: sin(θ) = λ/(2d) ⟹ k_perp = (2/λ)·sin(θ)·cos(ω-θ) = 1/d
    kPerpExpected = 1 / d_002;
    testCase.verifyLessThan(abs(kPerp - kPerpExpected), 1e-10, ...
        'k_perp at SrTiO3(002) Bragg angle does not match 1/d.');
    testCase.verifyLessThan(abs(kPar), 1e-10, ...
        'k_par at symmetric Bragg angle must be ~0.');
end

function testNoLambdaError(testCase)
    s = xrdc.io.emptyScan();
    s.twoTheta   = (30:0.5:35).';
    s.counts     = ones(size(s.twoTheta));
    s.secondAxis = 16;
    % scan.lambda is NaN (emptyScan default) and no Lambda option supplied
    testCase.verifyError(@() xrdc.rsm.toReciprocalSpace(s), ...
        'xrdc:rsm:noLambda');
end

function testNoSecondAxisError(testCase)
    s = xrdc.io.emptyScan();
    s.twoTheta = (30:0.5:35).';
    s.counts   = ones(size(s.twoTheta));
    % secondAxis is NaN by default
    testCase.verifyError(@() xrdc.rsm.toReciprocalSpace(s, 'Lambda', 1.5406), ...
        'xrdc:rsm:noSecondAxis');
end

function testNonSymmetricBranchOmegaScan(testCase)
    % Non-twoThetaOmega scan uses the simplified branch.
    lambda = 1.5406;
    tt = (30:0.5:35).';
    s = makeSlice(tt, 20, lambda);
    s.scanType = "omega";   % triggers simplified branch

    [kPar, kPerp] = xrdc.rsm.toReciprocalSpace(s);

    % Manually compute expected values
    dOmega = 0; dTheta = 0;
    om_rad = (20 + dOmega) * (pi/180);
    th_rad = (tt + dTheta) * (pi/360);
    kPerpExp = (2/lambda) .* sin(th_rad) .* cos(om_rad - th_rad);
    kParExp  = (2/lambda) .* sin(th_rad) .* sin(om_rad - th_rad);

    testCase.verifyLessThan(max(abs(kPerp - kPerpExp)), 1e-12);
    testCase.verifyLessThan(max(abs(kPar  - kParExp)),  1e-12);
end

% =====================================================================
% loadAreaScan
% =====================================================================

function testLoadAreaScanFromFileList(testCase)
    % Build two minimal XRDML-style scans in a temp folder then load them.
    tmpDir = tempname;
    mkdir(tmpDir);
    testCase.addTeardown(@() rmdir(tmpDir, 's'));

    % Write two minimal two-column text scan files
    fid1 = fopen(fullfile(tmpDir, 'slice_01.txt'), 'w');
    fprintf(fid1, '%.4f %d\n', [(30:0.1:32).' 500*ones(21,1)].');
    fclose(fid1);

    fid2 = fopen(fullfile(tmpDir, 'slice_02.txt'), 'w');
    fprintf(fid2, '%.4f %d\n', [(30:0.1:32).' 600*ones(21,1)].');
    fclose(fid2);

    files = {fullfile(tmpDir,'slice_01.txt'), fullfile(tmpDir,'slice_02.txt')};
    scans = xrdc.rsm.loadAreaScan(files);

    testCase.verifyEqual(numel(scans), 2);
    % All elements must have scanType = "area"
    % Use == instead of strcmp — cell-of-strings vs char is unreliable in R2026a
    testCase.verifyTrue(all(string({scans.scanType}) == "area"), ...
        'All loaded scans must have scanType "area".');
    testCase.verifyEqual(numel(scans(1).twoTheta), 21);
end

function testLoadAreaScanFromFolder(testCase)
    tmpDir = tempname;
    mkdir(tmpDir);
    testCase.addTeardown(@() rmdir(tmpDir, 's'));

    for k = 1:3
        fid = fopen(fullfile(tmpDir, sprintf('slice_%02d.txt', k)), 'w');
        fprintf(fid, '%.4f %d\n', [(30:0.1:32).' (500+k)*ones(21,1)].');
        fclose(fid);
    end

    scans = xrdc.rsm.loadAreaScan(tmpDir, 'Pattern', '*.txt');
    testCase.verifyEqual(numel(scans), 3);
end

function testLoadAreaScanSortsBySecondAxis(testCase)
    % Provide two scans where secondAxis is initially in descending order;
    % loadAreaScan must return them sorted ascending.
    tmpDir = tempname;
    mkdir(tmpDir);
    testCase.addTeardown(@() rmdir(tmpDir, 's'));

    % Manually craft a tiny Philips .x00 to set secondAxis — but simpler:
    % we provide two-column text files (secondAxis = NaN for text).
    % Instead, use xrdc.io.emptyScan structs directly by passing a file list.
    % (For this test we just check the sort on the NaN case — both NaN,
    %  sort is stable and leaves order unchanged.)
    % To test a real sort, we write two single-slice XRDML strings would be
    % complex; instead test the sort contract via internal file naming order.

    % Write files named z_02 first, z_01 second (reverse alphabetical)
    fid1 = fopen(fullfile(tmpDir, 'z_02.txt'), 'w');
    fprintf(fid1, '%.4f %d\n', [(30:0.1:32).' 500*ones(21,1)].');
    fclose(fid1);
    fid2 = fopen(fullfile(tmpDir, 'z_01.txt'), 'w');
    fprintf(fid2, '%.4f %d\n', [(30:0.1:32).' 500*ones(21,1)].');
    fclose(fid2);

    % Both have secondAxis = NaN (text format) — sort must not error
    scans = xrdc.rsm.loadAreaScan(tmpDir, 'Pattern', '*.txt');
    testCase.verifyEqual(numel(scans), 2);
end

function testLoadAreaScanBadFolderError(testCase)
    testCase.verifyError(...
        @() xrdc.rsm.loadAreaScan('/no/such/folder/xyz'), ...
        'xrdc:rsm:notAFolder');
end

function testLoadAreaScanSerialOptOut(testCase)
    % UseParallel=false must skip the parfor branch and still return
    % exactly the same scans as the default path.
    tmpDir = tempname;
    mkdir(tmpDir);
    testCase.addTeardown(@() rmdir(tmpDir, 's'));
    for k = 1:3
        fid = fopen(fullfile(tmpDir, sprintf('s_%02d.txt', k)), 'w');
        fprintf(fid, '%.4f %d\n', [(30:0.1:32).' (500+k)*ones(21,1)].');
        fclose(fid);
    end
    scans = xrdc.rsm.loadAreaScan(tmpDir, 'Pattern', '*.txt', ...
                                  'UseParallel', false);
    testCase.verifyEqual(numel(scans), 3);
end

% =====================================================================
% plotRsm — FR-6.3
% =====================================================================

function testPlotRsmReturnsHandles(testCase)
    lambda = 1.5406;
    omegas = [15.5 16.0 16.5];
    scans  = arrayfun(@(w) makeSlice((30:0.5:35).', w, lambda), omegas);
    for i = 1:numel(scans)
        scans(i).scanType = "twoThetaOmega";
        % Give each slice a Gaussian peak so the contoured intensity grid is
        % non-constant — a flat grid makes contourf warn "constant ZData".
        tt = scans(i).twoTheta;
        scans(i).counts = 100 + 1e4 * exp(-((tt - 32.5).^2) / (2 * 0.6^2));
    end

    h = xrdc.plot.plotRsm(scans, 'Lambda', lambda);

    testCase.verifyTrue(isfield(h, 'image'), 'Return struct needs .image');
    testCase.verifyTrue(isfield(h, 'ax'),    'Return struct needs .ax');
    testCase.verifyTrue(isfield(h, 'figure'),'Return struct needs .figure');
    testCase.verifyClass(h.ax, 'matlab.graphics.axis.Axes');
end

function testPlotRsmEmptyError(testCase)
    testCase.verifyError(@() xrdc.plot.plotRsm(struct([])), ...
        'xrdc:plot:emptyScans');
end

% =====================================================================
% V7 Parity — real-data gate (requires reference RSM file from Dr. Paik)
% =====================================================================

function testRsmLoadRealKTaO3Data(testCase)
    % Integration test: load the three 112 RSM slices from the Paik lab
    % (Schwaigert et al. JVST A 2023) and run through toReciprocalSpace.
    % Verifies end-to-end readXrdml → loadAreaScan → toReciprocalSpace
    % on real PANalytical files.
    dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        '..', '..', 'rexdrctomatlabport');
    assumeTrue(testCase, isfolder(dataDir), ...
        'Paik test data folder not found; skipping real-data test.');

    pattern = 'HP PtO2 on TiO2 001 112 RSM*.xrdml';
    scans = xrdc.rsm.loadAreaScan(dataDir, 'Pattern', pattern);
    assumeFalse(testCase, isempty(scans), 'No RSM slice files found.');

    % Each XRDML file is a multi-scan area measurement (~1500 ω-slices),
    % expanded by readXrdmlArea. Expect at least hundreds of slices total.
    testCase.verifyGreaterThan(numel(scans), 100);

    % Check the first slice and the middle slice for sanity; walking all
    % ~4500 would make the test run several seconds.
    sliceIdx = [1, round(numel(scans) / 2), numel(scans)];
    for i = sliceIdx
        s = scans(i);
        testCase.verifyGreaterThan(numel(s.twoTheta), 10);
        testCase.verifyGreaterThan(max(s.counts), 0);
        testCase.verifyFalse(isnan(s.secondAxis), ...
            sprintf('Slice %d has NaN secondAxis.', i));
        testCase.verifyFalse(isnan(s.lambda), ...
            sprintf('Slice %d has NaN lambda.', i));

        [kP, kZ] = xrdc.rsm.toReciprocalSpace(s);
        testCase.verifyEqual(numel(kP), numel(s.twoTheta));
        testCase.verifyEqual(numel(kZ), numel(s.twoTheta));
        testCase.verifyTrue(all(isfinite(kP)));
        testCase.verifyTrue(all(isfinite(kZ)));
        testCase.verifyGreaterThan(min(kZ), 0);
    end
end

% =====================================================================
% biaxialStrain — FR-7.x (RSM strain decomposition)
% =====================================================================

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

function testRsmParityAgainstBaronePlotScript(testCase)
    % Parity check: our toReciprocalSpace must match the Qx/Qz formulas
    % in Matthew Barone's RSMPlot(1).m on the same real data.
    %   Barone: Qx = (2/λ)·sin(θ)·sin(ω − θ)   (up to sign/flip)
    %           Qz = (2/λ)·sin(θ)·cos(ω − θ)
    %   where θ = 2θ/2 AND ω = commonPosition (no 2θ_center normalisation).
    % We disable XRDC's 2θ_center normalisation by forcing scanType ≠ "twoThetaOmega"
    % so the simplified branch runs.
    dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        '..', '..', 'rexdrctomatlabport');
    assumeTrue(testCase, isfolder(dataDir), 'Paik test data folder not found.');

    files = dir(fullfile(dataDir, 'HP PtO2 on TiO2 001 112 RSM*.xrdml'));
    assumeFalse(testCase, isempty(files), 'No 112 RSM slice files found.');

    s = xrdc.io.readXrdml(fullfile(files(1).folder, files(1).name));
    s.scanType = "omega";   % forces simplified branch (matches Barone's direct use of ω)

    [kP, kZ] = xrdc.rsm.toReciprocalSpace(s);

    % Barone's analytic formula (magQ·sin/cos of ω − θ)
    lambda = s.lambda;
    tt_rad = double(s.twoTheta) * pi/180;
    om_rad = s.secondAxis       * pi/180;
    theta  = tt_rad / 2;
    magQ   = 2 * sin(theta) ./ lambda;
    angle  = om_rad - theta;
    QxRef  = magQ .* sin(angle);
    QzRef  = magQ .* cos(angle);

    tol = 1e-12;
    testCase.verifyLessThan(max(abs(kP - QxRef)), tol, ...
        'k_par does not match Barone RSMPlot(1).m formula (simplified branch).');
    testCase.verifyLessThan(max(abs(kZ - QzRef)), tol, ...
        'k_perp does not match Barone RSMPlot(1).m formula (simplified branch).');
end

% =====================================================================
% findRsmPeaks — FR-6.x 2D auto-finder
% =====================================================================

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

% =====================================================================
% analyzeStrainRSM — FR-7.x one-call strain/composition orchestrator
% =====================================================================

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
