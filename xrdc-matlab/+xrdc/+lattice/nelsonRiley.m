function result = nelsonRiley(twoTheta, latticeValues)
%NELSONRILEY  Nelson–Riley extrapolation of a lattice parameter.
%   result = xrdc.lattice.nelsonRiley(twoTheta, latticeValues)
%
%   Linearly regresses a_i vs. the Nelson–Riley function
%
%        NR(θ) = cos²θ / sin θ  +  cos²θ / θ_deg
%
%   and reports the intercept (the extrapolated lattice parameter at θ = 90°,
%   i.e. where NR → 0). The mixed radians/degrees in the two terms is a
%   historical convention — do not "fix" it.
%
%   Inputs
%     twoTheta       : vector of 2θ in degrees for the selected peaks
%     latticeValues  : vector of per-peak lattice parameter estimates in Å,
%                      already corrected for order (e.g. a_i = d_i · m_i).
%
%   Output (struct)
%     .a0            : intercept = extrapolated lattice parameter (Å)
%     .slope         : slope of the NR regression (Å per NR-unit)
%     .a0SE          : standard error of the intercept (Å, closed-form)
%     .slopeSE       : standard error of the slope
%     .nrX           : NR function values for each peak
%     .residuals     : latticeValues − (a0 + slope·nrX)
%     .rSquared      : coefficient of determination
%
%   Reference: xrdc3.pas line 284+, CalculateNelsonRiley. See also the
%   regression uncertainty formula in ALGORITHM_SPEC.md §6.3.

    arguments
        twoTheta       (:,1) double {mustBePositive}
        latticeValues  (:,1) double {mustBePositive}
    end

    if numel(twoTheta) ~= numel(latticeValues)
        error('xrdc:lattice:sizeMismatch', ...
            'twoTheta and latticeValues must have the same length.');
    end

    n = numel(twoTheta);
    if n < 2
        error('xrdc:lattice:tooFewPoints', ...
            'Nelson–Riley needs at least 2 peaks (got %d).', n);
    end

    theta    = deg2rad(twoTheta / 2);   % θ in radians
    thetaDeg = twoTheta / 2;            % θ in degrees (Delphi's convention)

    % Nelson–Riley x-axis
    nrX = cos(theta).^2 ./ sin(theta) + cos(theta).^2 ./ thetaDeg;

    % closed-form OLS (matches Geradenanpassung in xrdc1.pas)
    Sx  = sum(nrX);
    Sy  = sum(latticeValues);
    Sxx = sum(nrX.^2);
    Sxy = sum(nrX .* latticeValues);

    denom = n * Sxx - Sx^2;
    slope = (n * Sxy - Sx * Sy) / denom;
    a0    = (Sxx * Sy - Sx * Sxy) / denom;

    % residuals + R²
    fitted    = a0 + slope * nrX;
    residuals = latticeValues - fitted;
    ssRes = sum(residuals.^2);
    ssTot = sum((latticeValues - mean(latticeValues)).^2);
    rSquared = 1 - ssRes / max(ssTot, eps);

    % Closed-form OLS standard errors (ALGORITHM_SPEC.md §6.3).
    %
    % Note on Delphi parity: xrdc3.pas:281 writes
    %     deltay0 = sqrt((s2x*s)/(n*(n-2)*(n*s2x-sqr(sx))))
    % which is off by a factor of 1/n from the textbook OLS intercept variance
    %     Var(b0) = σ² · Σxᵢ² / (n · Sxx)    with σ² = RSS/(n-2)
    % i.e. the Delphi expression underestimates the SE by √n. We use the
    % correct textbook formula here (CLAUDE.md "match the algorithm, not the
    % implementation") and do *not* reproduce the Delphi scaling. The slope
    % SE has no Delphi counterpart.
    if n >= 3
        rss     = ssRes;
        a0SE    = sqrt((Sxx * rss) / ((n - 2) * denom));
        slopeSE = sqrt((n   * rss) / ((n - 2) * denom));
    else
        a0SE    = NaN;
        slopeSE = NaN;
    end

    result = struct( ...
        'a0',        a0, ...
        'slope',     slope, ...
        'a0SE',      a0SE, ...
        'slopeSE',   slopeSE, ...
        'nrX',       nrX, ...
        'residuals', residuals, ...
        'rSquared',  rSquared);
end
