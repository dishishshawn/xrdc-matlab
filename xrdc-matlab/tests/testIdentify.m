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
