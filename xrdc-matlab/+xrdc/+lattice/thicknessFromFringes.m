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
%     .slopeSinThetaPerIndex : fitted slope (for debugging)
%     .residualsDeg     : residuals in 2θ (degrees) from the fit

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
    p = polyfit(idx, sinTh, 1);
    m = p(1);
    b = p(2);
    thicknessFitNm = lambdaNm / (2 * m);

    % uncertainty on slope → uncertainty on thickness (propagation)
    residSinTh = sinTh - (m * idx + b);
    if n >= 3
        s = sqrt(sum(residSinTh.^2) / (n - 2));      % residual std dev
        Sxx = sum((idx - mean(idx)).^2);
        slopeSE = s / sqrt(Sxx);
        thicknessFitSeNm = (lambdaNm / 2) * slopeSE / m^2;
    else
        thicknessFitSeNm = NaN;
    end

    % convert residuals back to 2θ (degrees) for human inspection
    predictedTwoTheta = 2 * rad2deg(asin(m * idx + b));
    residualsDeg = twoThetaFringes - predictedTwoTheta;

    result = struct( ...
        'thicknessNm',           thicknessNm, ...
        'thicknessFitNm',        thicknessFitNm, ...
        'thicknessFitSeNm',      thicknessFitSeNm, ...
        'slopeSinThetaPerIndex', m, ...
        'residualsDeg',          residualsDeg);
end
