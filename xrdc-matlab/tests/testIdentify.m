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
