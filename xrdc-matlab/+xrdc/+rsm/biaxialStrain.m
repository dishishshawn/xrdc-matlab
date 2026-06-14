function [a0, epsPar, epsPerp] = biaxialStrain(aPar, aPerp, factor)
%BIAXIALSTRAIN  Relaxed lattice parameter and biaxial strain of a (001) film.
%   [a0, epsPar, epsPerp] = xrdc.rsm.biaxialStrain(aPar, aPerp, factor)
%
%   Closed-form biaxial-strain decomposition for a (001)-oriented film
%   (Hooke's law, sigma_zz = 0 at the free surface). Given the measured
%   in-plane parameter aPar, out-of-plane parameter aPerp, and the elastic
%   factor f (= 2nu/(1-nu) or 2 c13/c33; see xrdc.lattice.elasticFactor):
%
%       a0      = (aPerp + f*aPar) / (1 + f)        relaxed (pseudocubic) param
%       epsPar  = (aPar  - a0) / a0                 in-plane strain  (<0 = compr.)
%       epsPerp = (aPerp - a0) / a0                 out-of-plane strain
%
%   These satisfy epsPerp/epsPar = -f and the reconstruction identity
%   a = a0*(1+eps) exactly. For an intrinsically tetragonal film (PTO/PZT)
%   a0 is the pseudocubic strain-model average, NOT a physical relaxed cubic
%   constant -- see docs/SCIENTIFIC_ASSUMPTIONS.md.
%
%   Inputs are scalar lattice parameters in Angstrom; factor is dimensionless.

    arguments
        aPar   (1,1) double {mustBePositive}
        aPerp  (1,1) double {mustBePositive}
        factor (1,1) double {mustBeNonnegative}
    end
    a0      = (aPerp + factor*aPar) / (1 + factor);
    epsPar  = (aPar  - a0) / a0;
    epsPerp = (aPerp - a0) / a0;
end
