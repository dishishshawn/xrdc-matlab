function energyEv = lambdaToEnergy(lambda)
%LAMBDATOENERGY  Photon energy in eV from wavelength in Å.
%   energyEv = xrdc.lattice.lambdaToEnergy(lambda)
%
%   See xrdc.lattice.energyToLambda for derivation.
%
%   Input
%     lambda    : wavelength in Å (scalar or array)
%
%   Output
%     energyEv  : photon energy in eV (same size as input)

    arguments
        lambda (:,:) double {mustBePositive}
    end

    h = 6.626068e-34;   % J·s
    c = 299792458;      % m/s
    e = 1.602e-19;      % C
    energyEv = 1e10 * (h * c) ./ (lambda * e);
end
