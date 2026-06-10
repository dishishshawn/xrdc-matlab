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

function testIdentifyUnsortedSubstrateCMeas(tc)
    % Substrate cMeas must pair each claimed d with ITS matched order,
    % not assume descending-d input (regression for order mispairing).
    lambda = 1.5406;
    ttSub  = xrdc.lattice.dToTwoTheta(3.905 ./ (1:4), lambda);
    tt     = ttSub([3 1 4 2]);                 % fixed shuffle
    R = xrdc.lattice.identifyMaterial(tt(:), lambda);
    tc.verifyTrue(R.substrate.found);
    tc.verifyEqual(R.substrate.cMeas, 3.905, 'AbsTol', 2e-3);
end

function testIdentifySubstrateOnlyNoSeries(tc)
    % Substrate-only scan: no film series, no crash on the empty table.
    lambda = 1.5406;
    ttSub  = xrdc.lattice.dToTwoTheta(3.905 ./ (1:4), lambda);
    R = xrdc.lattice.identifyMaterial(ttSub(:), lambda);
    tc.verifyTrue(R.substrate.found);
    tc.verifyEqual(height(R.series), 0);
    tc.verifyEmpty(R.unassigned);
end

function testIdentifyBadSubstrateErrors(tc)
    % PbTiO3 has role "film" - not usable as a declared substrate.
    tc.verifyError(@() xrdc.lattice.identifyMaterial([22.7; 46.5], 1.5406, ...
        Substrate="PbTiO3"), 'xrdc:lattice:badSubstrate');
end

function testIdentifyUnassignedStrayPeak(tc)
    % STO series + one stray peak at d = 2.31 A. As a single it is tried
    % at every order l = 1..4 -> c in {2.31, 4.62, 6.93, 9.24} A, none of
    % which lands in any candidate window (+- pad ~0.06 A):
    %   SrRuO3 hull [3.930, 3.951], PbTiO3 [4.1511, 4.152],
    %   PZT [4.146, 4.2613]. Must land in R.unassigned, not in a series.
    lambda  = 1.5406;
    ttSub   = xrdc.lattice.dToTwoTheta(3.905 ./ (1:4), lambda);
    ttStray = xrdc.lattice.dToTwoTheta(2.31, lambda);
    R = xrdc.lattice.identifyMaterial([ttSub(:); ttStray], lambda);
    tc.verifyTrue(R.substrate.found);
    tc.verifyEqual(height(R.series), 0);
    tc.verifyEqual(R.unassigned, ttStray, 'AbsTol', 1e-9);
end

function testIdentifyTieBreakParsimony(tc)
    % c = 4.1510, 1e-4 BELOW PTO's pseudomorphic prediction (4.15109):
    % PTO's raw score dips just under 1.0 while PZT's pseudomorphic
    % inversion keeps PZT at exactly 1.0. Quantized-score tie-break must
    % still rank fixed-stoichiometry PTO first (parsimony), with the
    % ambiguity flag set.
    tt = stoPlusFilmPeaks(4.1510, 1:3);
    R  = xrdc.lattice.identifyMaterial(tt, 1.5406);
    tc.verifyEqual(height(R.series), 1);
    tc.verifyEqual(R.series.bestMatch(1), "PbTiO3");
    cand = R.series.candidates{1};
    tc.verifyTrue(any(string({cand.name}) == "PZT"));
    tc.verifyTrue(any(R.series.flags{1} == "ambiguous"));
end
