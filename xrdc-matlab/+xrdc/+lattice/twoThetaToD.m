function d = twoThetaToD(twoTheta, lambda)
%TWOTHETATOD  d-spacing from 2θ via Bragg's law.
%   d = xrdc.lattice.twoThetaToD(twoTheta, lambda)
%
%   Bragg's law with n = 1:   d = λ / (2 sin θ)
%
%   Inputs
%     twoTheta : 2θ in degrees (scalar or array)
%     lambda   : wavelength in Å (scalar)
%
%   Output
%     d        : d-spacing in Å (same size as twoTheta)
%
%   See also xrdc.lattice.dToTwoTheta, xrdc.lattice.energyToLambda.

    arguments
        twoTheta  (:,:) double
        lambda    (1,1) double {mustBePositive}
    end

    theta = deg2rad(twoTheta) / 2;
    d = lambda ./ (2 * sin(theta));
end
