function twoTheta = dToTwoTheta(d, lambda)
%DTOTWOTHETA  2θ from d-spacing via Bragg's law.
%   twoTheta = xrdc.lattice.dToTwoTheta(d, lambda)
%
%   Bragg's law with n = 1:   2θ = 2·arcsin(λ / (2d))
%
%   Inputs
%     d       : d-spacing in Å (scalar or array)
%     lambda  : wavelength in Å (scalar)
%
%   Output
%     twoTheta : 2θ in degrees (same size as d); NaN where λ/(2d) > 1
%                (reflection impossible).
%
%   See also xrdc.lattice.twoThetaToD.

    arguments
        d       (:,:) double
        lambda  (1,1) double {mustBePositive}
    end

    ratio = lambda ./ (2 * d);
    twoTheta = nan(size(d));
    valid = ratio <= 1 & ratio >= -1;
    twoTheta(valid) = 2 * rad2deg(asin(ratio(valid)));
end
