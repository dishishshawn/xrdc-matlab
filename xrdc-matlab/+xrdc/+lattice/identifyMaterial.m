function R = identifyMaterial(peaks, lambda, opts)
%IDENTIFYMATERIAL  Identify material(s) from (00l) peak positions.
%   R = xrdc.lattice.identifyMaterial(peaks, lambda, Name=Value)
%
%   Pipeline (see docs/superpowers/specs/2026-06-09-material-id-design.md):
%   ghost filter -> claim declared substrate series -> harmonic grouping
%   -> per-series least-squares c -> strain-aware candidate matching.
%   Candidates are scored against the interval [c_bulk, c_pred] where
%   c_pred is the pseudomorphic-on-substrate prediction
%   (eps_perp = -f*eps_par, f = 2nu/(1-nu) or 2*c13/c33); a film sits
%   anywhere on that interval depending on relaxation. ALL accepted
%   candidates are returned ranked - never a silent single winner.
%
%   Inputs
%     peaks  : one of
%              - vector of peak 2-theta positions (degrees)
%              - struct array from xrdc.peaks.findPeaks (.twoTheta, .counts)
%              - scan struct (vector .twoTheta/.counts) -> findPeaks is
%                run with auto prominence and 0.2-deg separation
%     lambda : wavelength in A (default Cu Kalpha1 1.5406)
%
%   Name-Value
%     Substrate    (default "SrTiO3") - declared substrate (name or alias)
%     SubstrateTol (default 0.001)    - relative d window to claim substrate
%     Tolerance    (default 0.005)    - harmonic-series grouping tolerance
%     MaxOrder     (default 4)
%     WindowPad    (default 0.015)    - candidate window pad (relative c)
%     GhostFilter  (default true)
%
%   Output R (struct):
%     .substrate : .name, .found, .twoTheta, .cMeas, .cRef
%     .series    : table - cMeas, cSigma, orders, twoTheta, bestMatch,
%                  bestScore, candidates (cell of struct arrays with
%                  .name .score .misfit .cBulk .cPred .strainVsBulk
%                  .relaxation .x), flags (cell of string arrays)
%     .unassigned: 2-theta of peaks in no series and matching nothing
%     .ghosts    : table from xrdc.peaks.filterGhostPeaks
%     .notes     : string array (incl. the PZT strain-confounding caveat)
%
%   See also xrdc.lattice.loadMaterials, xrdc.lattice.groupHarmonicSeries,
%            xrdc.peaks.filterGhostPeaks, xrdc.peaks.findPeaks.

    arguments
        peaks
        lambda (1,1) double {mustBePositive} = 1.5406
        opts.Substrate    (1,1) string = "SrTiO3"
        opts.SubstrateTol (1,1) double {mustBePositive} = 0.001
        opts.Tolerance    (1,1) double {mustBePositive} = 0.005
        opts.MaxOrder     (1,1) double {mustBeInteger, mustBePositive} = 4
        opts.WindowPad    (1,1) double {mustBePositive} = 0.015
        opts.GhostFilter  (1,1) logical = true
    end

    PZT_CAVEAT = "Composition from (00l) alone is strain-confounded; the x" + ...
        " estimate assumes a fully pseudomorphic film. Deconvolving strain" + ...
        " and composition needs the in-plane parameter (RSM or asymmetric" + ...
        " reflection).";

    [tt, counts] = normalizeInput(peaks);
    if isempty(tt)
        error('xrdc:lattice:noPeaks', 'No peaks supplied or detected.');
    end

    % --- 1. ghost filter -------------------------------------------------
    ghosts = table();
    if opts.GhostFilter
        [keep, ghosts] = xrdc.peaks.filterGhostPeaks(tt, counts, lambda);
        tt = tt(keep);   % counts are not consumed past this point
    end

    % --- 2. substrate: declared input, confirmed not discovered ----------
    sub = xrdc.lattice.loadMaterials(opts.Substrate);
    if ~ismember(string(sub.role), ["substrate", "both"])
        error('xrdc:lattice:badSubstrate', ...
            '%s has role "%s" - not usable as a substrate.', sub.name, sub.role);
    end
    d = xrdc.lattice.twoThetaToD(tt, lambda);
    cSub = sub.c;
    isSubPeak = false(size(d));
    matchedL = zeros(size(d));   % order l claimed by each peak (0 = none)
    for l = 1:opts.MaxOrder
        dPred = cSub / l;
        [err, j] = min(abs(d - dPred) / dPred);
        if err <= opts.SubstrateTol && ~isSubPeak(j)
            isSubPeak(j) = true;
            matchedL(j) = l;
        end
    end
    % Pair each claimed d with ITS matched order (indexing by the same
    % mask) - pairing by loop order is wrong for unsorted input.
    R.substrate = struct('name', string(sub.name), 'found', any(isSubPeak), ...
        'twoTheta', tt(isSubPeak), 'cRef', cSub, ...
        'cMeas', meanOrNaN(d(isSubPeak) .* matchedL(isSubPeak)));
    if ~R.substrate.found
        warning('xrdc:lattice:substrateNotFound', ...
            ['Declared substrate %s: no (00l) series found within %.2f%%. ', ...
             'Proceeding with all peaks unassigned.'], ...
            sub.name, 100 * opts.SubstrateTol);
    end

    % --- 3. group remaining peaks into harmonic series --------------------
    ttF = tt(~isSubPeak);  dF = d(~isSubPeak);
    [S, singleIdx] = xrdc.lattice.groupHarmonicSeries(dF, ...
        Tolerance=opts.Tolerance, MaxOrder=opts.MaxOrder);

    % --- 4. strain-aware naming, with order-ambiguity hypotheses ----------
    M = xrdc.lattice.loadMaterials();
    films = M(ismember(string({M.role}), ["film", "both"]));
    rows = cell(0, 8);
    unassigned = [];
    for s = S
        hyps = orderHypotheses(s, opts.MaxOrder);
        bestH = []; bestCand = []; bestScoreH = -Inf; flags = strings(0,1);
        anyAcceptedCount = 0;
        for h = hyps
            cand = rankCandidates(films, h.c, sub.a, opts.WindowPad);
            if ~isempty(cand), anyAcceptedCount = anyAcceptedCount + 1; end
            hScore = -Inf;
            if ~isempty(cand), hScore = cand(1).score; end
            if hScore > bestScoreH || isempty(bestH)
                bestH = h; bestCand = cand; bestScoreH = hScore;
            end
        end
        if bestH.kind ~= "asIs", flags(end+1) = bestH.kind; end %#ok<AGROW>
        if anyAcceptedCount > 1, flags(end+1) = "orderAmbiguous"; end %#ok<AGROW>
        if numel(bestCand) >= 2 && bestCand(2).score >= bestCand(1).score - 0.2
            flags(end+1) = "ambiguous"; %#ok<AGROW>
        end
        bestName = ""; bestScore = NaN;
        if ~isempty(bestCand)
            bestName = string(bestCand(1).name); bestScore = bestCand(1).score;
        end
        rows(end+1, :) = {bestH.c, s.cSigma * bestH.scale, bestH.orders, ...
            ttF(s.members), bestName, bestScore, {bestCand}, {flags}}; %#ok<AGROW>
    end

    % Singles: try every order against the database; best (l, candidate) wins.
    for j = singleIdx(:).'
        bestCand = []; bestC = NaN; bestL = NaN;
        for l = 1:opts.MaxOrder
            cand = rankCandidates(films, dF(j) * l, sub.a, opts.WindowPad);
            if ~isempty(cand) && (isempty(bestCand) || cand(1).score > bestCand(1).score)
                bestCand = cand; bestC = dF(j) * l; bestL = l;
            end
        end
        if isempty(bestCand)
            unassigned(end+1) = ttF(j); %#ok<AGROW>
        else
            flags = "singlePeak";
            if numel(bestCand) >= 2 && bestCand(2).score >= bestCand(1).score - 0.2
                flags(end+1) = "ambiguous"; %#ok<AGROW>
            end
            rows(end+1, :) = {bestC, NaN, bestL, ttF(j), ...
                string(bestCand(1).name), bestCand(1).score, {bestCand}, {flags}}; %#ok<AGROW>
        end
    end

    seriesVars = {'cMeas','cSigma','orders','twoTheta','bestMatch', ...
        'bestScore','candidates','flags'};
    if isempty(rows)
        % cell2table on empty rows would leave 'candidates' non-cell and
        % break downstream cellfun calls; build the empty table explicitly.
        R.series = cell2table(cell(0, 8), 'VariableNames', seriesVars);
    else
        R.series = cell2table(rows, 'VariableNames', seriesVars);
    end
    R.unassigned = unassigned(:);
    R.ghosts = ghosts;
    R.lambda = lambda;
    R.notes = strings(0, 1);
    if height(R.series) > 0 && ...
            any(cellfun(@(c) any(string({c.name}) == "PZT"), R.series.candidates))
        R.notes(end+1) = PZT_CAVEAT;
    end
end

% ====================== local functions ======================

function [tt, counts] = normalizeInput(peaks)
    if isnumeric(peaks)
        tt = peaks(:); counts = nan(size(tt));
    elseif isstruct(peaks) && (numel(peaks) > 1 || ...
           (isscalar(peaks) && isscalar(peaks.twoTheta)))
        tt = [peaks.twoTheta].'; counts = [peaks.counts].';
    elseif isstruct(peaks) && isscalar(peaks)   % scan struct
        % Auto log-domain prominence; 0.2° separation merges the Kα1/Kα2
        % substrate split before harmonic grouping.
        pk = xrdc.peaks.findPeaks(peaks, 'MinSeparation', 0.2);
        if isempty(pk), tt = []; counts = []; return, end
        tt = [pk.twoTheta].'; counts = [pk.counts].';
    else
        tt = []; counts = [];
    end
end

function m = meanOrNaN(v)
    if isempty(v), m = NaN; else, m = mean(v); end
end

function hyps = orderHypotheses(s, maxOrder)
%ORDERHYPOTHESES  The series as grouped, doubled, and (if even-only) halved.
    hyps = struct('kind', "asIs", 'c', s.c, 'orders', s.orders, 'scale', 1);
    if 2 * max(s.orders) <= maxOrder
        hyps(end+1) = struct('kind', "orderDoubled", 'c', 2 * s.c, ...
            'orders', 2 * s.orders, 'scale', 2);
    end
    if s.evenOnly
        hyps(end+1) = struct('kind', "orderHalved", 'c', s.c / 2, ...
            'orders', s.orders / 2, 'scale', 0.5);
    end
end

function cand = rankCandidates(films, cMeas, aSub, pad)
%RANKCANDIDATES  Accepted candidates for a measured c, best first.
    cand = struct('name', {}, 'score', {}, 'misfit', {}, 'cBulk', {}, ...
        'cPred', {}, 'strainVsBulk', {}, 'relaxation', {}, 'x', {}, ...
        'hasComposition', {});
    for e = films
        [hull, cBulk, cPred, x] = candidateWindow(e, cMeas, aSub);
        if cMeas >= hull(1) && cMeas <= hull(2)
            misfit = 0;
        else
            misfit = min(abs(cMeas - hull)) / e.c;
        end
        if misfit > pad, continue, end
        denom = cBulk - cPred;
        if abs(denom) < 1e-4, relax = NaN;
        else, relax = (cMeas - cPred) / denom; end
        cand(end+1) = struct('name', string(e.name), ...
            'score', max(0, 1 - misfit / pad), 'misfit', misfit, ...
            'cBulk', cBulk, 'cPred', cPred, ...
            'strainVsBulk', (cMeas - cBulk) / cBulk, ...
            'relaxation', relax, 'x', x, ...
            'hasComposition', ~isempty(e.composition)); %#ok<AGROW>
    end
    if isempty(cand), return, end
    % Rank: quantized score desc, fixed-composition first (parsimony),
    % then |cMeas-cPred| asc (quantized at 1e-3 A).
    % The score key is quantized at 0.05: a composition-model entry's
    % pseudomorphic inversion gives score 1.0 EXACTLY across its whole
    % hull (cPred = cMeas by construction), so any sub-mA noise that
    % nudges a fixed entry's score below 1.0 would otherwise hand the
    % win to the composition model on the raw-score key. Scores within
    % 0.05 are treated as tied and parsimony (fixed stoichiometry)
    % decides (spec: dilute PZT can imitate PTO arbitrarily well, so
    % PTO must win exact AND near ties). Raw scores are still reported,
    % and the 0.2-window "ambiguous" flag is unchanged.
    key = [-round([cand.score] / 0.05); double([cand.hasComposition]); ...
           round(abs(cMeas - [cand.cPred]) / 1e-3)].';
    [~, order] = sortrows(key);
    cand = cand(order);
end

function [hull, cBulk, cPred, x] = candidateWindow(e, cMeas, aSub)
%CANDIDATEWINDOW  [min,max] of {c_bulk, c_pred} (over composition for PZT).
    f = strainFactor(e);
    x = NaN;
    if isempty(e.composition)
        cBulk = e.c;
        epsPar = (aSub - e.a) / e.a;
        cPred  = e.c * (1 - f * epsPar);
        hull   = [min(cBulk, cPred), max(cBulk, cPred)];
    else
        xg = linspace(min(e.composition.x), max(e.composition.x), 105);
        ab = interp1(e.composition.x, e.composition.a, xg);
        cb = interp1(e.composition.x, e.composition.c, xg);
        cp = cb .* (1 - f * (aSub - ab) ./ ab);
        hull = [min([cb cp]), max([cb cp])];
        % pseudomorphic inversion: cp is monotonic in x (a(x) increasing)
        if cMeas >= min(cp) && cMeas <= max(cp)
            x = interp1(cp, xg, cMeas);
            cBulk = interp1(xg, cb, x);
            cPred = cMeas;                       % on the pseudomorphic line
        else
            cBulk = e.c;                          % nominal entry
            epsPar = (aSub - e.a) / e.a;
            cPred  = e.c * (1 - f * epsPar);
        end
    end
end

function f = strainFactor(e)
    f = xrdc.lattice.elasticFactor(e);
end
