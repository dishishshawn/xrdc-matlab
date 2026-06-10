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
    [keep, ~] = xrdc.peaks.filterGhostPeaks(tt, nan(2,1), 1.5406);
    tc.verifyEqual(keep, [true; true]);   % nothing removed blind
end

function testGhostFilterWLineAndMixedNaN(tc)
    % --- W-Lalpha1 ghost removal ---
    % STO (002): d = 3.905/2; W-Lalpha1 ghost of that reflection.
    lambda  = 1.5406;
    d002    = 3.905/2;
    parent2t = xrdc.lattice.dToTwoTheta(d002, lambda);
    ghost2t  = xrdc.lattice.dToTwoTheta(d002, 1.4763);   % W-Lalpha1
    tt_w     = [ghost2t; parent2t];
    counts_w = [5e4; 1e6];   % ghost is 5% of parent — well below MaxRatio
    [keep_w, ghosts_w] = xrdc.peaks.filterGhostPeaks(tt_w, counts_w, lambda);
    tc.verifyEqual(keep_w, [false; true], ...
        'W-Lalpha1 ghost should be removed');
    tc.verifyEqual(height(ghosts_w), 1);
    tc.verifyEqual(ghosts_w.ghostLambda(1), 1.4763, 'AbsTol', 1e-9);

    % --- Mixed NaN/real: NaN-count peak must never be flagged ---
    % Setup: parent at STO (002), Kbeta ghost position occupied by a
    % NaN-count peak (unmeasured).  Even though it sits at the ghost
    % position it must survive because its intensity cannot be ranked.
    kbeta2t  = xrdc.lattice.dToTwoTheta(d002, 1.3922);
    tt_mix   = [kbeta2t; parent2t];          % NaN-peak at ghost pos, parent
    counts_mix = [NaN; 1e6];
    % No warning should be emitted (not all-NaN).
    tc.verifyWarningFree( ...
        @() xrdc.peaks.filterGhostPeaks(tt_mix, counts_mix, lambda));
    [keep_mix, ghosts_mix] = xrdc.peaks.filterGhostPeaks(tt_mix, counts_mix, lambda);
    tc.verifyTrue(keep_mix(1), ...
        'NaN-count peak at ghost position must not be flagged');
    tc.verifyTrue(keep_mix(2), 'Parent must remain');
    tc.verifyEqual(height(ghosts_mix), 0, ...
        'No ghosts should be recorded when the candidate has NaN counts');
end

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
