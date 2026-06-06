function result = thicknessQuadratic(twoThetaFringes, lambda)
%THICKNESSQUADRATIC  Modified-Bragg (Kiessig) film thickness from fringes.
%   result = xrdc.xrr.thicknessQuadratic(twoThetaFringes, lambda)
%
%   The full Kiessig relation between fringe positions and film thickness is
%
%        sin²(θ_n) = sin²(θ_c) + (n · λ / 2t)²
%
%   not the linear approximation sin(θ_n) ≈ sin(θ_c) + n · λ/(2t) that
%   xrdc.lattice.thicknessFromFringes uses. The linear form is only valid
%   when sin(θ_c) << n · λ/(2t), which breaks down for thin films (small
%   t → wide fringe spacing → comparable to sin(θ_c)). For typical
%   30-50 nm PLD oxide films on Cu Kα radiation, the linear formula
%   over-estimates t by ~5%.
%
%   With N detected fringes at indices i = 1..N and the true (but
%   unknown) first-fringe index n_0, the data follow
%
%        sin²(θ_i) = α + β·i + γ·i²
%
%   with γ = (λ/2t)², β = 2·γ·n_0, α = sin²(θ_c) + γ·n_0². A polynomial
%   fit therefore yields thickness, fringe-index offset, AND critical
%   edge in one shot, with no need to assume θ_c separately.
%
%   Inputs
%     twoThetaFringes : 2θ positions of detected fringes (degrees, N ≥ 3)
%     lambda          : wavelength in Å
%
%   Output (struct)
%     .thicknessNm        Film thickness (nm)
%     .thicknessSeNm      1-σ SE on thickness
%     .thicknessCi95Nm    [lo hi] 95% CI on thickness
%     .twoThetaC          Recovered critical edge (2θ degrees)
%     .firstFringeIndex   Implied true index of the first detected fringe
%                          (small float; non-integer = our index labelling
%                          is offset by this much from the physical count)
%     .rSquared           Goodness of fit
%     .residualsDeg       Residuals back in 2θ (degrees) for inspection
%     .fitBackend         "cfit" when Curve Fitting Toolbox available;
%                          otherwise "polyfit".

    arguments
        twoThetaFringes (:,1) double
        lambda          (1,1) double {mustBePositive}
    end

    twoThetaFringes = sort(twoThetaFringes);
    n = numel(twoThetaFringes);
    if n < 3
        error('xrdc:xrr:tooFewFringes', ...
            'Quadratic Kiessig fit needs ≥3 fringes, got %d.', n);
    end

    sinTh    = sin(twoThetaFringes * pi / 360);    % sin(θ); 2θ → θ via /2
    y        = sinTh.^2;
    idx      = (1:n).';
    lambdaNm = lambda / 10;

    useCfit = ~isempty(which('fit'));
    if useCfit
        try
            [cf, gof] = fit(idx, y, 'poly2');
            g = cf.p1;  b = cf.p2;  a = cf.p3;     % p1=γ p2=β p3=α
            rSquared   = gof.rsquare;
            fitBackend = "cfit";

            ci = confint(cf, 0.95);
            gLo = ci(1, 1); gHi = ci(2, 1);
            % t = λ/(2·√γ): non-monotonic transform handled by bracketing.
            tLo = lambdaNm / (2 * sqrt(max(gHi, eps)));
            tHi = lambdaNm / (2 * sqrt(max(gLo, eps)));
            thicknessCi95Nm = sort([tLo, tHi]);
            slopeSE  = (gHi - gLo) / (2 * 1.96);
            thicknessSeNm = abs(lambdaNm / 4) * slopeSE / max(g, eps)^(3/2);
        catch ME
            warning('xrdc:xrr:cfitFailed', ...
                'Curve Fitting fit() failed (%s). Falling back to polyfit.', ME.message);
            useCfit = false;
        end
    end

    if ~useCfit
        p = polyfit(idx, y, 2);
        g = p(1); b = p(2); a = p(3);
        yHat = polyval(p, idx);
        resid = y - yHat;
        if n >= 4
            sigma2 = sum(resid.^2) / (n - 3);
            V = [idx.^2, idx, ones(n,1)];
            covP = sigma2 * inv(V.' * V);     %#ok<MINV> — small system
            seG = sqrt(max(covP(1,1), 0));
            thicknessSeNm = abs(lambdaNm / 4) * seG / max(g, eps)^(3/2);
        else
            thicknessSeNm = NaN;
        end
        thicknessCi95Nm = [NaN NaN];
        rSquared = 1 - sum(resid.^2) / max(sum((y - mean(y)).^2), eps);
        fitBackend = "polyfit";
    end

    if g <= 0
        error('xrdc:xrr:nonPhysicalFit', ...
            'Quadratic Kiessig fit returned γ ≤ 0 — not a fringe pattern.');
    end

    thicknessNm    = lambdaNm / (2 * sqrt(g));
    nZero          = b / (2 * g);                  % first-fringe offset
    sinSqThetaC    = max(a - g * nZero^2, 0);
    twoThetaC      = 2 * asin(sqrt(sinSqThetaC)) * 180 / pi;

    yHat = a + b * idx + g * idx.^2;
    predictedTwoTheta = 2 * asin(sqrt(max(yHat, 0))) * 180 / pi;
    residualsDeg = twoThetaFringes - predictedTwoTheta;

    result = struct( ...
        'thicknessNm',       thicknessNm, ...
        'thicknessSeNm',     thicknessSeNm, ...
        'thicknessCi95Nm',   thicknessCi95Nm, ...
        'twoThetaC',         twoThetaC, ...
        'firstFringeIndex',  nZero, ...
        'rSquared',          rSquared, ...
        'residualsDeg',      residualsDeg, ...
        'fitBackend',        fitBackend);
end
