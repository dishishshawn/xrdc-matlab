function result = thicknessFromFringes(twoThetaFringes, lambda)
%THICKNESSFROMFRINGES  Film thickness from Kiessig / Laue fringe positions.
%   result = xrdc.lattice.thicknessFromFringes(twoThetaFringes, lambda)
%
%   For N fringes at 2θ₁ < 2θ₂ < … < 2θ_N, the film thickness is
%
%        t = (N − 1) · λ / (2 · (sin θ_N − sin θ_1))
%
%   The full-precision version fits sin(θ_i) vs. fringe index i with a
%   straight line; the slope is λ / (2 t) which gives a thickness
%   uncertainty for free. Both are returned.
%
%   Matches xrdc15.pas (EstimateThickness).
%
%   Inputs
%     twoThetaFringes : vector of 2θ fringe positions in degrees (≥ 2 entries)
%     lambda          : wavelength in Å (scalar)
%
%   Output (struct)
%     .thicknessNm      : thickness from N-fringe formula (nm)
%     .thicknessFitNm   : thickness from linear fit of sin θ vs. index (nm)
%     .thicknessFitSeNm : 1-sigma uncertainty on the fit thickness (nm)
%     .thicknessCi95Nm  : [lo, hi] 95% CI on the fit thickness (nm)
%                         (Curve Fitting Toolbox path; NaN otherwise)
%     .slopeSinThetaPerIndex : fitted slope (for debugging)
%     .rSquared         : fit goodness (1.0 = perfect linear sin θ vs. n)
%     .residualsDeg     : residuals in 2θ (degrees) from the fit
%     .fitBackend       : "cfit" if the Curve Fitting Toolbox was used,
%                         "polyfit" otherwise.
%
%   The Curve Fitting Toolbox is used for confidence-interval reporting
%   when available; otherwise the routine falls back to polyfit and the
%   1-sigma SE is reported (ci field stays NaN).

    arguments
        twoThetaFringes (:,1) double
        lambda          (1,1) double {mustBePositive}
    end

    twoThetaFringes = sort(twoThetaFringes);
    n = numel(twoThetaFringes);
    if n < 2
        error('xrdc:lattice:tooFewFringes', ...
            'Need at least 2 fringe positions, got %d.', n);
    end

    theta = deg2rad(twoThetaFringes / 2);
    sinTh = sin(theta);

    lambdaNm = lambda / 10;         % Å → nm
    % N-fringe formula (original Delphi)
    thicknessNm = (n - 1) * lambdaNm / (2 * (sinTh(end) - sinTh(1)));

    % Linear fit: sin θ_i = m · i + b, where m = λ / (2t)
    idx = (1:n).';

    useCfit = ~isempty(which('fit')) && exist('fitoptions', 'file') == 2;
    if useCfit && n >= 2
        % Curve Fitting Toolbox path: gives a real confidence interval and
        % an R² alongside the slope. With only 2 fringes the CI collapses
        % to NaN (interval undefined) but the slope is still valid.
        try
            [cf, gof] = fit(idx, sinTh, 'poly1');
            m = cf.p1;
            b = cf.p2;
            rSquared = gof.rsquare;
            fitBackend = "cfit";

            ci = confint(cf, 0.95);   % 2×2: rows = [lo;hi], cols = [p1,p2]
            mLo = ci(1, 1);
            mHi = ci(2, 1);
            thicknessCi95Nm = sort([lambdaNm / (2 * mLo), lambdaNm / (2 * mHi)]);
            % 1-σ SE from CI half-width / 1.96 (Gaussian approx)
            slopeSE = (mHi - mLo) / (2 * 1.96);
            thicknessFitSeNm = (lambdaNm / 2) * slopeSE / m^2;
        catch ME
            warning('xrdc:lattice:cfitFailed', ...
                'Curve Fitting fit() failed (%s). Falling back to polyfit.', ME.message);
            useCfit = false;
        end
    end

    if ~useCfit
        % polyfit fallback (no Curve Fitting Toolbox).
        p = polyfit(idx, sinTh, 1);
        m = p(1);
        b = p(2);
        residSinTh = sinTh - (m * idx + b);
        if n >= 3
            s = sqrt(sum(residSinTh.^2) / (n - 2));
            Sxx = sum((idx - mean(idx)).^2);
            slopeSE = s / sqrt(Sxx);
            thicknessFitSeNm = (lambdaNm / 2) * slopeSE / m^2;
            rSquared = 1 - sum(residSinTh.^2) / max(sum((sinTh - mean(sinTh)).^2), eps);
        else
            thicknessFitSeNm = NaN;
            rSquared = NaN;
        end
        thicknessCi95Nm = [NaN, NaN];
        fitBackend = "polyfit";
    end

    thicknessFitNm = lambdaNm / (2 * m);

    % convert residuals back to 2θ (degrees) for human inspection
    predictedTwoTheta = 2 * rad2deg(asin(m * idx + b));
    residualsDeg = twoThetaFringes - predictedTwoTheta;

    result = struct( ...
        'thicknessNm',           thicknessNm, ...
        'thicknessFitNm',        thicknessFitNm, ...
        'thicknessFitSeNm',      thicknessFitSeNm, ...
        'thicknessCi95Nm',       thicknessCi95Nm, ...
        'slopeSinThetaPerIndex', m, ...
        'rSquared',              rSquared, ...
        'residualsDeg',          residualsDeg, ...
        'fitBackend',            fitBackend);
end
