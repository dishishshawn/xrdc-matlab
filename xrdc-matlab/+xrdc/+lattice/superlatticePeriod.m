function result = superlatticePeriod(twoThetaSatellites, lambda)
%SUPERLATTICEPERIOD  Superlattice period Λ from satellite-peak positions.
%   result = xrdc.lattice.superlatticePeriod(twoThetaSatellites, lambda)
%
%   A superlattice [A/B]×N produces satellite reflections evenly spaced in
%   sin θ around the average Bragg peak. Adjacent satellite orders obey the
%   same finite-period relation as Kiessig / Laue fringes:
%
%        Λ = (M − 1) · λ / (2 · (sin θ_M − sin θ_1))
%
%   so the bilayer period is recovered with the identical math used for film
%   thickness. This wraps xrdc.lattice.thicknessFromFringes and relabels its
%   output as a period (the linear sin θ vs. satellite-order fit gives the
%   period uncertainty for free).
%
%   This returns the SUPERLATTICE PERIOD (the A+B bilayer repeat), NOT the
%   individual A and B layer thicknesses — separating those requires full
%   dynamical / optical (Parratt) modelling, which is out of scope here.
%   Total stack thickness comes from XRR (xrdc.xrr.analyzeFringes).
%
%   Inputs
%     twoThetaSatellites : vector of satellite 2θ positions (degrees, ≥ 2),
%                          adjacent orders around a single main peak.
%     lambda             : wavelength in Å (scalar).
%
%   Output (struct)
%     .periodNm        Period Λ from the M-satellite endpoint formula (nm)
%     .periodFitNm     Period from the sin θ vs. order linear fit (nm)
%     .periodFitSeNm   1-σ uncertainty on the fit period (nm)
%     .periodCi95Nm    [lo, hi] 95% CI on the fit period (nm); NaN without
%                      the Curve Fitting Toolbox
%     .rSquared        Linear-fit goodness (1.0 = perfectly periodic)
%     .nSatellites     Number of satellite positions used
%     .residualsDeg    Per-satellite residuals in 2θ (degrees)
%     .fitBackend      "cfit" or "polyfit"
%
%   See also XRDC.LATTICE.THICKNESSFROMFRINGES, XRDC.XRR.ANALYZEFRINGES.

    arguments
        twoThetaSatellites (:,1) double
        lambda             (1,1) double {mustBePositive}
    end

    if numel(twoThetaSatellites) < 2
        error('xrdc:lattice:tooFewSatellites', ...
            'Need at least 2 satellite positions, got %d.', ...
            numel(twoThetaSatellites));
    end

    % The satellite spacing in sin θ is identical to fringe spacing, so the
    % thickness routine recovers the period directly — reuse it (DRY) and
    % relabel the fields as a period.
    t = xrdc.lattice.thicknessFromFringes(twoThetaSatellites, lambda);

    result = struct( ...
        'periodNm',      t.thicknessNm, ...
        'periodFitNm',   t.thicknessFitNm, ...
        'periodFitSeNm', t.thicknessFitSeNm, ...
        'periodCi95Nm',  t.thicknessCi95Nm, ...
        'rSquared',      t.rSquared, ...
        'nSatellites',   numel(twoThetaSatellites), ...
        'residualsDeg',  t.residualsDeg, ...
        'fitBackend',    t.fitBackend);
end
