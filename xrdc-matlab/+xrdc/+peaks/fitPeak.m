function result = fitPeak(scan, window, options)
%FITPEAK  Nonlinear fit of a single diffraction peak.
%   result = xrdc.peaks.fitPeak(scan, [tthMin, tthMax])
%   result = xrdc.peaks.fitPeak(scan, [tthMin, tthMax], Name, Value, ...)
%
%   Modern replacement for the brute-force grid search in
%   TForm8.DoFit (xrdc9.pas:161). Uses `lsqcurvefit` by default; a
%   `'Method','bruteforce'` option keeps the old algorithm available for
%   reproducing legacy analyses. See ALGORITHM_SPEC §5.
%
%   Input
%     scan   : scan struct with .twoTheta, .counts.
%     window : [twoThetaMin, twoThetaMax] — the fit window around a
%              single peak. Must contain exactly one peak.
%
%   Name/Value options
%     'Shape'          "lorentz" (default) | "gauss" | "pseudoVoigt"
%     'Background'     "linear"  (default) | "none"
%                      "linear" fits y = bga·(x − tthMin) + bgb in
%                      addition to the peak shape.
%     'BgFraction'     (default 0.15)
%                      Fraction of the window on each side used to
%                      initialise the linear background (Delphi
%                      bgposleft / bgposright, xrdc9.pas:188).
%     'InitialGuess'   struct with optional fields
%                       .center, .fwhm, .amplitude, .eta
%                      Any missing field is auto-initialised from the
%                      data.
%     'Method'         "lsqcurvefit" (default) | "bruteforce"
%                      bruteforce replicates xrdc9.pas:161+.
%     'MaxIter'        (default 400) — passed to lsqcurvefit.
%     'TolFun'         (default 1e-10)
%
%   Returns (struct)
%     .shape           shape used
%     .twoTheta        fitted peak centre (degrees)
%     .fwhm            fitted FWHM (degrees)
%     .amplitude       peak height above background
%     .eta             Lorentz fraction (pseudo-Voigt only; NaN else)
%     .bga, .bgb       linear background slope/intercept (0 if 'none')
%     .paramSE         struct with SE for each fitted parameter
%                      (Jacobian-derived; NaN under 'bruteforce')
%     .rss             residual sum of squares
%     .rmse            residual RMS (in counts)
%     .rSquared        coefficient of determination
%     .xFit, .yFit     2θ and model evaluated on the fit window
%     .residuals       counts - yFit at the data points
%     .dataX, .dataY   the points used in the fit
%     .method          method used
%
%   Requires the Optimization Toolbox when 'Method' = 'lsqcurvefit'.

    arguments
        scan                          (1,1) struct
        window                        (1,2) double
        options.Shape                 (1,1) string {mustBeMember( ...
            options.Shape, ["lorentz","gauss","pseudoVoigt"])} = "lorentz"
        options.Background            (1,1) string {mustBeMember( ...
            options.Background, ["linear","none"])} = "linear"
        options.BgFraction            (1,1) double {mustBePositive, ...
            mustBeLessThan(options.BgFraction, 0.5)} = 0.15
        options.InitialGuess          (1,1) struct = struct()
        options.Method                (1,1) string {mustBeMember( ...
            options.Method, ["lsqcurvefit","bruteforce"])} = "lsqcurvefit"
        options.MaxIter               (1,1) double {mustBeInteger, mustBePositive} = 400
        options.TolFun                (1,1) double {mustBePositive} = 1e-10
    end

    if window(1) >= window(2)
        error('xrdc:peaks:badWindow', ...
            'window must be [min, max] with min < max.');
    end

    x = double(scan.twoTheta(:));
    y = double(scan.counts(:));
    mask = x >= window(1) & x <= window(2);
    if nnz(mask) < 5
        error('xrdc:peaks:tooFewPoints', ...
            'Fit window contains only %d points (need >= 5).', nnz(mask));
    end
    xi = x(mask);
    yi = y(mask);

    % Initial linear background (matches xrdc9.pas:188-209). Use the
    % first/last BgFraction of the window; take the mean count there.
    nw = numel(xi);
    nb = max(1, round(options.BgFraction * nw));
    yL = mean(yi(1:nb));
    yR = mean(yi(end-nb+1:end));
    xL = mean(xi(1:nb));
    xR = mean(xi(end-nb+1:end));
    if abs(xR - xL) < eps
        bga0 = 0; bgb0 = (yL + yR) / 2;
    else
        bga0 = (yR - yL) / (xR - xL);
        bgb0 = yL - bga0 * xL;
    end

    % Background-subtracted residual, peak height from its max
    yBg0   = bga0 * xi + bgb0;
    yPk    = yi - yBg0;
    [ampG, iMax] = max(yPk);
    centerG = xi(iMax);

    % FWHM guess: width where yPk drops to half of its max (or ~1/5 of
    % window if we can't find it).
    halfLevel = ampG / 2;
    above = yPk > halfLevel;
    if any(above)
        iAbove = find(above);
        fwhmG  = xi(iAbove(end)) - xi(iAbove(1));
        if fwhmG <= 0, fwhmG = (xi(end) - xi(1)) / 5; end
    else
        fwhmG  = (xi(end) - xi(1)) / 5;
    end

    guess = options.InitialGuess;
    if ~isfield(guess, 'center'),    guess.center    = centerG; end
    if ~isfield(guess, 'fwhm'),      guess.fwhm      = fwhmG;   end
    if ~isfield(guess, 'amplitude'), guess.amplitude = max(ampG, eps); end
    if ~isfield(guess, 'eta'),       guess.eta       = 0.5;     end

    % ---- Dispatch ----
    if options.Method == "lsqcurvefit"
        % Try lsqcurvefit; fall back to bruteforce if Optimization Toolbox missing
        try
            result = fitLsq(xi, yi, guess, options, bga0, bgb0);
        catch ME
            if contains(ME.message, 'Optimization Toolbox')
                warning('xrdc:peaks:noOptToolbox', ...
                    'Optimization Toolbox not available. Falling back to bruteforce method.');
                result = fitBruteforce(xi, yi, guess, options, bga0, bgb0);
            else
                rethrow(ME);
            end
        end
    else
        result = fitBruteforce(xi, yi, guess, options, bga0, bgb0);
    end

    result.dataX      = xi;
    result.dataY      = yi;
    result.xFit       = xi;
    result.yFit       = evalModel( ...
        result.twoTheta, result.fwhm, result.amplitude, ...
        result.eta, result.bga, result.bgb, options.Shape, xi);
    result.residuals  = yi - result.yFit;
    result.rss        = sum(result.residuals.^2);
    result.rmse       = sqrt(result.rss / numel(yi));
    result.rSquared   = 1 - result.rss / max(sum((yi - mean(yi)).^2), eps);
    result.shape      = options.Shape;
    result.method     = options.Method;
end

% =========================================================================

function r = fitLsq(xi, yi, guess, options, bga0, bgb0)
    shape = options.Shape;
    useBg = options.Background == "linear";
    [params0, lb, ub, meta] = buildParamVector(guess, bga0, bgb0, shape, useBg, xi);

    modelFn = @(p, x) evalPacked(p, x, shape, meta);

    opts = optimoptions('lsqcurvefit', ...
        'Display', 'off', ...
        'MaxIterations', options.MaxIter, ...
        'FunctionTolerance', options.TolFun, ...
        'StepTolerance', 1e-12);

    [p, resnorm, residuals, ~, ~, ~, J] = ...
        lsqcurvefit(modelFn, params0, xi, yi, lb, ub, opts); %#ok<ASGLU>

    [center, fwhm, amp, eta, bga, bgb] = unpack(p, meta, bga0, bgb0);

    % Jacobian-based parameter standard errors:
    %     σ² = RSS / (n - k)
    %     Cov = σ² · (J'·J)^-1
    n    = numel(yi);
    k    = numel(p);
    if n > k
        sigma2 = resnorm / (n - k);
        % J may be sparse; Matlab returns full J anyway for small problems
        JtJ = full(J' * J);
        covP = sigma2 * pinv(JtJ);
        seP  = sqrt(max(diag(covP), 0));
    else
        seP  = nan(size(p));
    end

    r = struct();
    r.twoTheta  = center;
    r.fwhm      = fwhm;
    r.amplitude = amp;
    r.eta       = eta;
    r.bga       = bga;
    r.bgb       = bgb;
    r.paramSE   = unpackSE(seP, meta);
end

function r = fitBruteforce(xi, yi, guess, options, bga0, bgb0)
    % Faithful port of xrdc9.pas:161+. Only supports Lorentz and Gauss
    % (no pseudo-Voigt in the Delphi original).
    if options.Shape == "pseudoVoigt"
        error('xrdc:peaks:bruteforceShape', ...
            'pseudoVoigt shape is not available under Method=bruteforce.');
    end

    bga = bga0; bgb = bgb0;
    dataResid = yi - (bga * xi + bgb);

    % Search grids — match xrdc9.pas:225-252.
    step   = median(diff(xi));
    xscaleMin = 1 / (xi(end) - xi(1));
    xscaleMax = 10 / max(step, eps);
    pkMax  = max(dataResid);
    yscaleMin = max(1, (pkMax - (mean(dataResid(1:max(1,round(0.1*numel(xi)))))) ) / 3);
    yscaleMax = max(1, pkMax * 2);
    x0Min  = xi(1);
    x0Max  = xi(end);

    devBest = Inf;
    devMin  = Inf;
    x0Fit = guess.center;
    xsFit = fwhmToXscale(guess.fwhm, options.Shape);
    ysFit = guess.amplitude;

    % Converge in passes, shrinking the search grid each time (20 steps
    % per parameter per pass, matching Delphi).
    while true
        devBest = devMin;
        devMin  = Inf;

        % exponential grids for the scales, linear for x0
        if xscaleMin <= 0, xscaleMin = eps; end
        xscaleStep = (xscaleMax / xscaleMin) ^ (1/20);
        yscaleStep = (yscaleMax / yscaleMin) ^ (1/20);
        x0Step     = (x0Max - x0Min) / 20;

        xs = xscaleMin;
        while xs <= xscaleMax
            ys = yscaleMin;
            while ys <= yscaleMax
                xx = x0Min;
                while xx <= x0Max
                    switch options.Shape
                        case "lorentz"
                            lf = ys ./ (1 + (xs * (xi - xx)).^2);
                        case "gauss"
                            lf = ys .* exp(-(xs * (xi - xx)).^2);
                    end
                    dev = sum((lf - dataResid).^2);
                    if dev < devMin
                        devMin = dev;
                        x0Fit  = xx;
                        xsFit  = xs;
                        ysFit  = ys;
                    end
                    xx = xx + x0Step;
                end
                ys = ys * yscaleStep;
            end
            xs = xs * xscaleStep;
        end

        % Shrink ranges around the best point (xrdc9.pas:300-305)
        xscaleMin = max(xscaleMin, xsFit / (xscaleStep^2));
        xscaleMax = min(xscaleMax, xsFit * (xscaleStep^2));
        yscaleMin = max(yscaleMin, ysFit / (yscaleStep^2));
        yscaleMax = min(yscaleMax, ysFit * (yscaleStep^2));
        x0Min     = max(x0Min, x0Fit - 2*x0Step);
        x0Max     = min(x0Max, x0Fit + 2*x0Step);

        if abs(devMin - devBest) < devBest / 1e10 && x0Step < 0.001
            break
        end
    end

    fwhmFit = xscaleToFwhm(xsFit, options.Shape);

    r = struct();
    r.twoTheta  = x0Fit;
    r.fwhm      = fwhmFit;
    r.amplitude = ysFit;
    r.eta       = NaN;
    r.bga       = bga;
    r.bgb       = bgb;
    r.paramSE   = struct('center', NaN, 'fwhm', NaN, 'amplitude', NaN, ...
                         'eta', NaN, 'bga', NaN, 'bgb', NaN);
end

% =========================================================================
% Shape evaluation (normalised to unit amplitude at the centre)
% =========================================================================

function y = evalShape(x, x0, fwhm, shape)
    fwhm = max(fwhm, eps);
    switch shape
        case "lorentz"
            g = fwhm / 2;
            y = g.^2 ./ ((x - x0).^2 + g.^2);
        case "gauss"
            sigma = fwhm / (2*sqrt(2*log(2)));
            y = exp(-(x - x0).^2 ./ (2 * sigma.^2));
        otherwise
            % pseudoVoigt handled separately
            y = zeros(size(x));
    end
end

function y = evalModel(x0, fwhm, amp, eta, bga, bgb, shape, x)
    switch shape
        case "lorentz"
            ypk = amp .* evalShape(x, x0, fwhm, "lorentz");
        case "gauss"
            ypk = amp .* evalShape(x, x0, fwhm, "gauss");
        case "pseudoVoigt"
            yL = evalShape(x, x0, fwhm, "lorentz");
            yG = evalShape(x, x0, fwhm, "gauss");
            eta = min(max(eta, 0), 1);
            ypk = amp .* (eta .* yL + (1 - eta) .* yG);
    end
    y = ypk + bga .* x + bgb;
end

% =========================================================================
% Parameter packing for lsqcurvefit
% =========================================================================

function [p, lb, ub, meta] = buildParamVector(guess, bga0, bgb0, shape, useBg, xi)
    % Parameter order: [center, fwhm, amplitude, (eta), (bga), (bgb)]
    p     = [guess.center, guess.fwhm, guess.amplitude];
    lb    = [xi(1),        eps,        0];
    ub    = [xi(end),      xi(end)-xi(1), Inf];
    names = {'center', 'fwhm', 'amplitude'};
    meta  = struct('shape', shape, 'hasEta', false, 'hasBg', useBg);
    if shape == "pseudoVoigt"
        p     = [p, guess.eta];
        lb    = [lb, 0];
        ub    = [ub, 1];
        names = [names, {'eta'}];
        meta.hasEta = true;
    end
    if useBg
        p     = [p, bga0, bgb0];
        lb    = [lb, -Inf, -Inf];
        ub    = [ub,  Inf,  Inf];
        names = [names, {'bga', 'bgb'}];
    end
    meta.names = names;
    p  = p(:).';
    lb = lb(:).';
    ub = ub(:).';
end

function y = evalPacked(p, x, shape, meta)
    center = p(1);
    fwhm   = p(2);
    amp    = p(3);
    idx    = 4;
    if meta.hasEta
        eta = p(idx); idx = idx + 1;
    else
        eta = NaN;
    end
    if meta.hasBg
        bga = p(idx); bgb = p(idx+1);
    else
        bga = 0; bgb = 0;
    end
    y = evalModel(center, fwhm, amp, eta, bga, bgb, shape, x);
end

function [center, fwhm, amp, eta, bga, bgb] = unpack(p, meta, bga0, bgb0)
    center = p(1);
    fwhm   = p(2);
    amp    = p(3);
    idx    = 4;
    if meta.hasEta
        eta = p(idx); idx = idx + 1;
    else
        eta = NaN;
    end
    if meta.hasBg
        bga = p(idx); bgb = p(idx+1);
    else
        bga = bga0; bgb = bgb0;
    end
end

function se = unpackSE(seP, meta)
    se = struct('center', NaN, 'fwhm', NaN, 'amplitude', NaN, ...
                'eta', NaN, 'bga', NaN, 'bgb', NaN);
    se.center    = seP(1);
    se.fwhm      = seP(2);
    se.amplitude = seP(3);
    idx = 4;
    if meta.hasEta
        se.eta = seP(idx); idx = idx + 1;
    end
    if meta.hasBg
        se.bga = seP(idx);
        se.bgb = seP(idx+1);
    end
end

% =========================================================================
% Conversions between Delphi xscale and modern FWHM
% =========================================================================

function xscale = fwhmToXscale(fwhm, shape)
    switch shape
        case "lorentz"
            xscale = 2 / max(fwhm, eps);
        case "gauss"
            xscale = 2 * sqrt(-log(0.5)) / max(fwhm, eps);
        otherwise
            xscale = 2 / max(fwhm, eps);
    end
end

function fwhm = xscaleToFwhm(xscale, shape)
    switch shape
        case "lorentz"
            fwhm = 2 / max(xscale, eps);
        case "gauss"
            fwhm = 2 * sqrt(-log(0.5)) / max(xscale, eps);
        otherwise
            fwhm = 2 / max(xscale, eps);
    end
end
