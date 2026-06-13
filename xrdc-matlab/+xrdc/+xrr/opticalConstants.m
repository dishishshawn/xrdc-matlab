function [delta, beta] = opticalConstants(material, density, lambda)
%OPTICALCONSTANTS  X-ray refractive-index decrement for a material at Cu Ka.
%   [delta, beta] = xrdc.xrr.opticalConstants(material, density, lambda)
%
%   Returns the dispersion (delta) and absorption (beta) of the complex
%   refractive index n = 1 - delta + i*beta, computed from a chemical
%   formula and mass density via the atomic-scattering sum
%     delta = (re*lambda^2/2pi) * sum_j N_j f1_j
%     beta  = (re*lambda^2/2pi) * sum_j N_j f2_j
%   with N_j the number density of element j (atoms/Angstrom^3).
%
%   Inputs
%     material : a materials.json name/alias (its .formula is used) OR a
%                chemical-formula string directly (e.g. "SrTiO3").
%     density  : mass density in g/cm^3.
%     lambda   : wavelength in Angstroms. v1 supports Cu Ka only.
%
%   Errors: xrdc:xrr:unsupportedEnergy, xrdc:xrr:unknownElement.
%   See also xrdc.xrr.reflectivityModel, xrdc.xrr.fitReflectivity.

    arguments
        material (1,1) string
        density  (1,1) double {mustBePositive}
        lambda   (1,1) double {mustBePositive}
    end
    CU_KA = 1.5406;  TOL = 0.02;
    if abs(lambda - CU_KA) > TOL
        error('xrdc:xrr:unsupportedEnergy', ...
            ['opticalConstants uses a Cu Ka (%.4f A) scattering table; ' ...
             'got lambda = %.4f A.'], CU_KA, lambda);
    end

    formula = resolveFormula(material);
    [elems, mult] = parseFormula(formula);
    data = loadScattering();

    re = 2.8179403262e-5;   % classical electron radius, Angstrom
    NA = 6.02214076e23;     % mol^-1

    molarMass = 0; f1sum = 0; f2sum = 0;
    for i = 1:numel(elems)
        el = char(elems(i));
        if ~isfield(data, el)
            error('xrdc:xrr:unknownElement', ...
                'No Cu Ka scattering data for element "%s".', el);
        end
        e = data.(el);
        molarMass = molarMass + mult(i) * e.mass;
        f1sum = f1sum + mult(i) * e.f1;
        f2sum = f2sum + mult(i) * e.f2;
    end

    % formula units per Angstrom^3:
    % rho[g/cm^3] * NA[mol^-1] / M[g/mol] = units/cm^3; /1e24 -> units/Angstrom^3
    nFU  = density * NA / molarMass / 1e24;
    pref = re * lambda^2 / (2*pi);
    delta = pref * nFU * f1sum;
    beta  = pref * nFU * f2sum;
end

% -------------------------------------------------------------------------

function formula = resolveFormula(material)
%RESOLVEFORMULA  A database name -> its .formula; otherwise treat as a formula.
    if ~isempty(regexp(char(material), '\d', 'once'))
        formula = material; return            % has a digit -> already a formula
    end
    try
        m = xrdc.lattice.loadMaterials(material);
        formula = string(m.formula);
    catch
        formula = material;                    % bare element symbol / formula
    end
end

function [elems, mult] = parseFormula(formula)
%PARSEFORMULA  "SrTiO3" -> ("Sr","Ti","O"), (1,1,3). Fractional counts OK.
    tok = regexp(char(formula), '([A-Z][a-z]?)(\d*\.?\d*)', 'tokens');
    n = numel(tok);
    elems = strings(n,1); mult = zeros(n,1);
    for i = 1:n
        elems(i) = string(tok{i}{1});
        c = tok{i}{2};
        if isempty(c), mult(i) = 1; else, mult(i) = str2double(c); end
    end
end

function data = loadScattering()
%LOADSCATTERING  Cached element f1/f2/mass table (Cu Ka).
    persistent cache
    if isempty(cache)
        p = fullfile(fileparts(mfilename('fullpath')), '..', '+data', ...
            'atomicScattering.json');
        cache = jsondecode(fileread(p)).elements;
    end
    data = cache;
end
