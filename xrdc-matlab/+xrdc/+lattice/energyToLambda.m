function lambda = energyToLambda(energyEv)
%ENERGYTOLAMBDA  X-ray wavelength in Å from photon energy in eV.
%   lambda = xrdc.lattice.energyToLambda(energyEv)
%
%   λ[Å] = 1e10 · (h·c) / (E[J]) where E[J] = E[eV] · e
%        = 1e10 · (h·c / e) / E[eV]
%        ≈ 12398.4197 / E[eV]        (Å / eV)
%
%   Matches the original XRDC Delphi constants (xrdc3.pas):
%     h = 6.626068e-34 J·s,  c = 299792458 m/s,  e = 1.602e-19 C.
%
%   Input
%     energyEv : photon energy in eV (scalar or array)
%
%   Output
%     lambda   : wavelength in Å (same size as input)
%
%   See also xrdc.lattice.lambdaToEnergy.

    arguments
        energyEv (:,:) double {mustBePositive}
    end

    h = 6.626068e-34;   % J·s
    c = 299792458;      % m/s
    e = 1.602e-19;      % C
    lambda = 1e10 * (h * c) ./ (energyEv * e);
end
