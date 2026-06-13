function result = fitReflectivity(scan, options)
%FITREFLECTIVITY  Fit a single-film XRR slab model (thickness/density/roughness).
%   result = xrdc.xrr.fitReflectivity(scan, Name=Value)
%
%   Fits film thickness, density, roughness, substrate roughness, an overall
%   scale, and a constant background to an XRR scan by minimising a log-space
%   residual of a Parratt model (xrdc.xrr.reflectivityModel). Seeded from
%   xrdc.xrr.analyzeFringes (thickness) and xrdc.xrr.findCriticalEdge
%   (density).
%
%   Optimiser. The log-reflectivity of a finite film is a fast-oscillating
%   function of thickness (every Kiessig fringe shifts as t changes), so a
%   finite-difference Levenberg–Marquardt step (lsqnonlin alone) stalls in a
%   spurious local minimum well above the true RSS — even when seeded at the
%   correct thickness. We therefore run a bounded, derivative-free
%   Nelder–Mead search (fminsearch) from a small spread of thickness seeds
%   and keep the lowest-RSS basin; this reliably finds the global minimum on
%   synthetic round-trips. When the Optimization Toolbox is present we then
%   polish that optimum with lsqnonlin to recover an analytic Jacobian for
%   the parameter standard errors (pinv-based, mirroring fitPeak — finite
%   even when the footprint/scale columns make J'J rank-deficient). Without
%   the toolbox the SEs come from a numeric finite-difference Jacobian.
%
%   Name/Value: Film ("SrTiO3"), Substrate ("SrTiO3"), UpperBound (5.0),
%   EdgeBuffer (0.05), Footprint (true), MaxIter (400).
%   Output fields: thicknessNm/SeNm, densityGcc/SeGcc, densityFraction,
%   filmRoughnessNm/SeNm, substrateRoughnessNm/SeNm, scale, background,
%   footprintDeg, rSquared, chiSq, rssSeed, rssFinal, modelCurve, residuals,
%   window, converged, method, seed.  Errors: xrdc:lattice:unknownMaterial.

    arguments
        scan                (1,1) struct
        options.Film        (1,1) string  = "SrTiO3"
        options.Substrate   (1,1) string  = "SrTiO3"
        options.UpperBound  (1,1) double  = 5.0
        options.EdgeBuffer  (1,1) double  = 0.05
        options.Footprint   (1,1) logical = true
        options.MaxIter     (1,1) double {mustBeInteger, mustBePositive} = 400
    end

    x = double(scan.twoTheta(:));
    y = double(scan.counts(:));
    lambda = scan.lambda;

    filmMat = xrdc.lattice.loadMaterials(options.Film);
    subMat  = xrdc.lattice.loadMaterials(options.Substrate);

    % --- Seeds ---------------------------------------------------------
    edge = xrdc.xrr.findCriticalEdge(x, y);
    fr   = xrdc.xrr.analyzeFringes(scan, UpperBound=options.UpperBound);
    tSeed = fr.thicknessNm; if ~isfinite(tSeed) || tSeed <= 0, tSeed = 20; end
    densSeed = densityFromEdge(edge.twoThetaC, filmMat, lambda);
    roughSeed = 0.3;  subRoughSeed = 0.3;

    win = [edge.twoThetaC + options.EdgeBuffer, options.UpperBound];
    inWin = x >= win(1) & x <= win(2) & y > 0;
    xw = x(inWin); yw = y(inWin);
    logY = log10(yw);

    scaleSeed = max(yw);
    bgSeed    = max(1, min(yw));
    fpSeed    = edge.twoThetaC;

    % p: 1 thickness(nm) 2 density(g/cc) 3 filmRough(nm) 4 subRough(nm)
    %    5 scale 6 background 7 footprintDeg
    p0 = [tSeed, densSeed, roughSeed, subRoughSeed, scaleSeed, bgSeed, fpSeed];
    lb = [1,    0.3*filmMat.densityBulk, 0,   0,   scaleSeed*1e-2, 0,   0];
    ub = [1e4,  1.3*filmMat.densityBulk, 10,  10,  scaleSeed*1e2,  max(yw), 2];
    if ~options.Footprint
        p0(7) = 0; lb(7) = 0; ub(7) = 0;
    end

    modelLog = @(p) log10(max( ...
        p(5) * xrdc.xrr.reflectivityModel(xw, ...
            buildLayers(p, filmMat, subMat), lambda, ...
            Footprint=options.Footprint, FootprintDeg=p(7)) + p(6), 1e-300));
    resid = @(p) modelLog(p) - logY;
    rssFun = @(p) sum(resid(clampParams(p, lb, ub)).^2);

    rssSeed = sum(resid(p0).^2);

    % --- Stage 1: multi-start bounded Nelder-Mead -----------------------
    % Thickness is the only badly multimodal coordinate, so we restart over
    % a spread of thickness seeds (the fringe seed plus a fixed fan) and keep
    % the lowest-RSS basin. Other coordinates start from their physical seeds.
    nmOpts = optimset('Display','off', ...
        'MaxIter', options.MaxIter*15, 'MaxFunEvals', options.MaxIter*15, ...
        'TolFun', 1e-8, 'TolX', 1e-8);
    tSeeds = unique(max(lb(1), [tSeed, 10, 15, 20, 25, 30, 35, 45, 60]));
    pBest = p0; rssBest = Inf;
    for ts = tSeeds
        pStart = p0; pStart(1) = ts;
        [pRaw, rss] = fminsearch(rssFun, pStart, nmOpts);
        pNM = clampParams(pRaw, lb, ub);
        if rss < rssBest
            rssBest = rss; pBest = pNM;
        end
    end

    % --- Stage 2: polish + Jacobian SE ---------------------------------
    useLsq = ~isempty(which('lsqnonlin'));
    if useLsq
        opts = optimoptions('lsqnonlin', 'Display','off', ...
            'MaxIterations', options.MaxIter, ...
            'MaxFunctionEvaluations', 200*numel(p0));
        [pPol, rssPol, ~, exitflag, ~, ~, J] = ...
            lsqnonlin(resid, pBest, lb, ub, opts);
        if rssPol <= rssBest          % keep the polish only if it improves
            pFit = pPol; rssFinal = rssPol;
        else
            pFit = pBest; rssFinal = rssBest;
            J = numericJacobian(resid, pFit);
        end
        converged = exitflag >= 0;     % NM already converged; lsqnonlin refines
        method = "neldermead+lsqnonlin";
        se = jacobianSE(J, resid(pFit));
    else
        pFit = pBest; rssFinal = rssBest;
        converged = true;
        method = "neldermead";
        se = jacobianSE(numericJacobian(resid, pFit), resid(pFit));
    end

    modelCurve = pFit(5) * xrdc.xrr.reflectivityModel(x, ...
        buildLayers(pFit, filmMat, subMat), lambda, ...
        Footprint=options.Footprint, FootprintDeg=pFit(7)) + pFit(6);
    residuals = log10(max(modelCurve(inWin),1e-300)) - logY;
    ssTot = sum((logY - mean(logY)).^2);
    rSq = 1 - sum(residuals.^2)/ssTot;

    result = struct( ...
        'thicknessNm', pFit(1), 'thicknessSeNm', se(1), ...
        'densityGcc',  pFit(2), 'densitySeGcc',  se(2), ...
        'densityFraction', pFit(2)/filmMat.densityBulk, ...
        'filmRoughnessNm', pFit(3), 'filmRoughnessSeNm', se(3), ...
        'substrateRoughnessNm', pFit(4), 'substrateRoughnessSeNm', se(4), ...
        'scale', pFit(5), 'background', pFit(6), 'footprintDeg', pFit(7), ...
        'rSquared', rSq, 'chiSq', rssFinal, ...
        'rssSeed', rssSeed, 'rssFinal', rssFinal, ...
        'modelCurve', modelCurve, 'residuals', residuals, ...
        'window', win, 'converged', converged, 'method', method, ...
        'seed', struct('thicknessNm',tSeed,'densityGcc',densSeed, ...
                       'roughnessNm',roughSeed));
end

% -------------------------------------------------------------------------

function layers = buildLayers(p, filmMat, subMat)
    layers = struct( ...
        'material',  {string(filmMat.name), string(subMat.name)}, ...
        'density',   {p(2), subMat.densityBulk}, ...
        'thickness', {p(1), Inf}, ...
        'roughness', {p(3), p(4)});
end

function dens = densityFromEdge(twoThetaC, filmMat, lambda)
%DENSITYFROMEDGE  Invert θc = sqrt(2δ) to a density seed (δ ∝ density).
    try
        thetaC = (twoThetaC/2) * pi/180;
        deltaObs = 0.5 * thetaC^2;
        [deltaBulk, ~] = xrdc.xrr.opticalConstants( ...
            string(filmMat.name), filmMat.densityBulk, lambda);
        dens = filmMat.densityBulk * (deltaObs / deltaBulk);
        if ~isfinite(dens) || dens <= 0, dens = filmMat.densityBulk; end
        dens = min(max(dens, 0.4*filmMat.densityBulk), 1.2*filmMat.densityBulk);
    catch
        dens = filmMat.densityBulk;
    end
end

function J = numericJacobian(resid, p)
%NUMERICJACOBIAN  Forward-difference Jacobian of the residual vector at p.
    r0 = resid(p);
    n  = numel(r0);
    k  = numel(p);
    J  = zeros(n, k);
    for j = 1:k
        h = max(1e-6, 1e-6 * abs(p(j)));
        pj = p; pj(j) = pj(j) + h;
        J(:,j) = (resid(pj) - r0) / h;
    end
end

function se = jacobianSE(J, r)
%JACOBIANSE  Parameter SEs from the Jacobian (mirrors fitPeak's pinv path).
%   Uses pinv so a rank-deficient J'J (footprint/scale columns are nearly
%   degenerate in-window) still yields finite SEs for the identifiable
%   parameters rather than NaN.
    J = full(J);
    dof = max(1, numel(r) - size(J,2));
    s2  = sum(r.^2) / dof;
    JtJ = J.'*J;
    C   = s2 * pinv(JtJ);
    se  = sqrt(max(diag(C), 0));
    se  = se(:).';
end

function pc = clampParams(p, lb, ub)
    pc = min(max(p, lb), ub);
end
