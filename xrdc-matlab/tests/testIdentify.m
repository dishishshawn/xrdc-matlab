function tests = testIdentify
%TESTIDENTIFY  Tests for material identification (materials.json,
%   filterGhostPeaks, groupHarmonicSeries, identifyMaterial).
    tests = functiontests(localfunctions);
end

% ---------- loadMaterials ----------

function testLoadMaterialsAll(tc)
    M = xrdc.lattice.loadMaterials();
    tc.verifyGreaterThanOrEqual(numel(M), 7);
    tc.verifyTrue(all(isfield(M, {'name','aliases','system','a','c','role','elastic','composition'})));
    names = string({M.name});
    tc.verifyTrue(all(ismember(["SrTiO3","SrRuO3","PbTiO3","PZT","LaAlO3","TiO2","WO3"], names)));
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

function testGhostFilterToleranceSparesNearbyRealPeak(tc)
    % Regression for the PositionTol default (0.15 -> 0.10 deg): a genuine
    % weak film peak 0.12 deg away from a predicted ghost position must
    % survive.  Synthetic boundary check (parent uses STO's d002 for
    % convenience); the real case it guards is S31, where PTO 002 sits
    % 0.121 deg from the W-Lalpha ghost position of SRO 002 and was
    % falsely removed at tol 0.15.
    lambda = 1.5406;
    d002   = 3.905/2;
    parent = xrdc.lattice.dToTwoTheta(d002, lambda);
    ghost  = xrdc.lattice.dToTwoTheta(d002, 1.4763);   % W-Lalpha1 position
    tt = [ghost + 0.12; parent];  counts = [5e4; 1e6];
    keep = xrdc.peaks.filterGhostPeaks(tt, counts, lambda);
    tc.verifyEqual(keep, [true; true], ...
        'Peak 0.12 deg from a ghost position must not be flagged');
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

function testLoadMaterialsHasXrrFields(tc)
    % XRR needs formula + bulk density on every material, and the new
    % substrate/film entries LAO, TiO2, WO3 must be present.
    M = xrdc.lattice.loadMaterials();
    tc.verifyTrue(all(isfield(M, {'formula','densityBulk'})), ...
        'formula/densityBulk missing from materials.json');
    names = string({M.name});
    for want = ["LaAlO3","TiO2","WO3"]
        tc.verifyTrue(any(names == want), sprintf('%s entry missing', want));
    end
    sto = xrdc.lattice.loadMaterials("SrTiO3");
    tc.verifyEqual(string(sto.formula), "SrTiO3");
    tc.verifyGreaterThan(sto.densityBulk, 4.5);   % STO ~5.12 g/cm^3
    tc.verifyLessThan(sto.densityBulk, 5.5);
end

% ---------- real-data validation (gated) ----------

function d = validationInputDir()
%Helper (not a test): validation data dir, anchored on this file's location
%   (the unittest runner cd's into tests/, so cwd-relative paths are fragile).
    d = fullfile(fileparts(mfilename('fullpath')), '..', '..', ...
        'validation', 'tushar', 'input');
end

function testIdentifyS25RealData(tc)
    p = fullfile(validationInputDir(), ...
        'TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_2theta_omega_05062026.txt');
    tc.assumeTrue(isfile(p), 'S25 validation scan not present');
    scan = xrdc.io.readScan(p);
    % Scan-struct input (auto findPeaks at 5% prominence) misses the weak
    % SRO film peak: SRO 002 is ~6e4 counts vs ~6.3e6 for STO 002 (~1% -
    % even the 1% threshold clips it by a hair).  Find peaks explicitly.
    pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts)*0.005);
    R = xrdc.lattice.identifyMaterial(pk, 1.5406, Substrate="SrTiO3");
    tc.verifyTrue(R.substrate.found);
    tc.verifyTrue(any(R.series.bestMatch == "SrRuO3"), ...
        'SRO film series not identified');
end

function testIdentifyS31Heterostructure(tc)
    p = fullfile(validationInputDir(), 'Heterostructure raw data', ...
        ['TR_S31_1_PTO on SRO_STO(100)_580c_150mT_and_700c_ 200mT', ...
         '_10500sh_5hz_2theta_omega_05202026.txt']);
    tc.assumeTrue(isfile(p), 'S31 validation scan not present');
    scan = xrdc.io.readScan(p);
    % The PTO film peaks are very weak relative to the substrate (PTO 002
    % is ~4.6e3 counts vs ~1.3e7 for STO 002, i.e. ~0.04%), so the default
    % 5% scan-struct prominence misses them.  Find peaks explicitly.
    pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts)*2e-4);
    R = xrdc.lattice.identifyMaterial(pk, 1.5406, Substrate="SrTiO3");
    tc.verifyTrue(R.substrate.found);
    % PTO must appear among candidates of some series (bestMatch may be
    % flagged ambiguous vs dilute PZT - that is correct behaviour).
    hasPTO = any(cellfun(@(c) ~isempty(c) && ...
        any(string({c.name}) == "PbTiO3"), R.series.candidates));
    tc.verifyTrue(hasPTO, 'PTO not among any series candidates');
end

function testElasticFactorNuAndC13C33(tc)
    % nu path: f = 2 nu / (1 - nu)
    sto = xrdc.lattice.loadMaterials("SrTiO3");   % nu = 0.232
    tc.verifyEqual(xrdc.lattice.elasticFactor(sto), ...
        2*0.232/(1-0.232), 'AbsTol', 1e-12);
    % c13/c33 path: f = 2 c13/c33
    e = struct('elastic', struct('c13', 70, 'c33', 90));
    tc.verifyEqual(xrdc.lattice.elasticFactor(e), 2*70/90, 'AbsTol', 1e-12);
    % error path: neither present
    bad = struct('elastic', struct());
    tc.verifyError(@() xrdc.lattice.elasticFactor(bad), 'xrdc:lattice:noElastic');
end

% ---------- KNOWN BUG: (00l)-only PTO-vs-PZT mis-ranking (gated) ----------
%
% On a real PbTiO3 film the measured out-of-plane c (~4.11 A) sits BELOW
% both PbTiO3's strain window ([cBulk 4.152, cPred 4.151]) and PZT's
% ([4.146, 4.261]) - the film has reduced tetragonality, c < bulk. The
% scorer ranks dilute PZT just above PbTiO3 purely because PZT's bulk c is
% 0.006 A lower (closer to 4.11), NOT because the film is Zr-doped. The
% series is correctly flagged 'ambiguous' and PbTiO3 is always the immediate
% runner-up, so the answer is hedged-wrong, never silently wrong.
%
% These are CHARACTERIZATION tests: they assert the current (wrong) ranking
% so it is visible in the suite. The physically-correct best match is
% PbTiO3. The real fix needs in-plane a from an asymmetric RSM (Feature A).
% WHEN ONE OF THESE FAILS, the scoring/RSM fix probably landed - confirm
% PbTiO3 now wins, then invert the bestMatch assertion (do not just delete).

function d = dataDumpDir()
%Helper (not a test): repo-root data dump, anchored on this file's location.
    d = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'data');
end

function s = ptoFilmSeries(tc, R)
%Helper (not a test): the series whose candidates include PbTiO3 (the film).
    idx = find(cellfun(@(c) ~isempty(c) && ...
        any(string({c.name}) == "PbTiO3"), R.series.candidates));
    tc.assertNotEmpty(idx, 'no series carries a PbTiO3 candidate');
    s = R.series(idx(1), :);
end

function testIdentifyS11KnownMisrankPtoVsPzt(tc)
    % S11: pure PbTiO3 film on STO. KNOWN-WRONG: PZT outranks PbTiO3.
    p = fullfile(dataDumpDir(), ...
        'TR_S11_PTO_STO(100)_580c_150mT_20000sh_5hz_2theta omega_04162026.txt');
    tc.assumeTrue(isfile(p), 'S11 data dump scan not present');
    scan = xrdc.io.readScan(p);
    pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts)*2e-4);
    R  = xrdc.lattice.identifyMaterial(pk, 1.5406, Substrate="SrTiO3");
    tc.verifyTrue(R.substrate.found);

    s = ptoFilmSeries(tc, R);
    tc.verifyGreaterThan(s.cMeas, 4.10);   % compressed-below-bulk film c
    tc.verifyLessThan(s.cMeas, 4.12);
    tc.verifyTrue(any(s.flags{1} == "ambiguous"), ...
        'film series should be flagged ambiguous');

    cand = s.candidates{1};
    names = string({cand.name});
    tc.verifyTrue(all(ismember(["PZT","PbTiO3"], names)), ...
        'both PZT and PbTiO3 must be candidates');
    % CHARACTERIZATION (known-wrong): PZT currently wins, PbTiO3 second.
    tc.verifyEqual(string(s.bestMatch), "PZT", ...
        'S11 ranking changed - PTO-vs-PZT fix may have landed; see header note');
    pzScore  = cand(names == "PZT").score;
    ptoScore = cand(names == "PbTiO3").score;
    tc.verifyGreaterThan(pzScore, ptoScore);
end

function testIdentifyS31KnownMisrankPtoFilm(tc)
    % S31 heterostructure: PTO film series mis-ranks the same way as S11,
    % while the SRO buffer + STO substrate series are identified correctly.
    p = fullfile(validationInputDir(), 'Heterostructure raw data', ...
        ['TR_S31_1_PTO on SRO_STO(100)_580c_150mT_and_700c_ 200mT', ...
         '_10500sh_5hz_2theta_omega_05202026.txt']);
    tc.assumeTrue(isfile(p), 'S31 validation scan not present');
    scan = xrdc.io.readScan(p);
    pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts)*2e-4);
    R  = xrdc.lattice.identifyMaterial(pk, 1.5406, Substrate="SrTiO3");
    tc.verifyTrue(R.substrate.found);
    % SRO buffer is correctly identified somewhere (sanity that the scan and
    % grouping are right, isolating the failure to the PTO film ranking).
    tc.verifyTrue(any(R.series.bestMatch == "SrRuO3"), ...
        'SRO buffer series not identified');

    s = ptoFilmSeries(tc, R);
    tc.verifyGreaterThan(s.cMeas, 4.09);
    tc.verifyLessThan(s.cMeas, 4.12);
    tc.verifyTrue(any(s.flags{1} == "ambiguous"), ...
        'film series should be flagged ambiguous');
    % CHARACTERIZATION (known-wrong): PZT currently wins, PbTiO3 second.
    tc.verifyEqual(string(s.bestMatch), "PZT", ...
        'S31 ranking changed - PTO-vs-PZT fix may have landed; see header note');
end
