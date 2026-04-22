function d = dSpacingFromHKL(h, k, l, lattice)
%DSPACINGFROMHKL  d-spacing for (h,k,l) reflection given crystal system.
%   d = xrdc.lattice.dSpacingFromHKL(h, k, l, lattice)
%
%   Inputs
%     h, k, l  : Miller indices (scalar or same-size arrays of integers)
%     lattice  : struct with fields (subset used depends on .system):
%                  .system  — 'cubic' | 'tetragonal' | 'orthorhombic' |
%                             'hexagonal' | 'rhombohedral' | 'monoclinic' |
%                             'triclinic'
%                  .a, .b, .c  (Å, required by the relevant system)
%                  .alpha, .beta, .gamma  (degrees, ditto)
%
%   Output
%     d        : d-spacing in Å, same size as h (NaN for forbidden (000)).
%
%   Formulas match xrdc1.pas:4341 (CalculateStructureLines).
%   Angles in the .alpha/.beta/.gamma fields are in degrees; they are
%   converted to radians internally. For cubic/tetragonal/orthorhombic/
%   hexagonal/monoclinic, missing angle fields default to the conventional
%   values (90° or 120° for hexagonal .gamma) so users can pass a
%   minimal struct.
%
%   See also xrdc.lattice.simulatePattern.

    arguments
        h         (:,:) double
        k         (:,:) double
        l         (:,:) double
        lattice   (1,1) struct
    end

    if ~isfield(lattice, 'system')
        error('xrdc:lattice:missingField', ...
            'lattice struct must have a .system field.');
    end

    % broadcast h, k, l to a common size
    [h, k, l] = broadcastTriplet(h, k, l);

    % normalise angle fields (degrees → radians), apply defaults
    defaultGamma = 90;
    if strcmpi(lattice.system, 'hexagonal')
        defaultGamma = 120;
    end
    alpha = getAngleRad(lattice, 'alpha', 90);
    beta  = getAngleRad(lattice, 'beta',  90);
    gamma = getAngleRad(lattice, 'gamma', defaultGamma); %#ok<NASGU>   (used in triclinic branch)

    % guard
    getA = @() mustHave(lattice, 'a');
    getB = @() mustHave(lattice, 'b');
    getC = @() mustHave(lattice, 'c');

    inv_d2 = nan(size(h));   % 1/d²

    switch lower(string(lattice.system))
        case "cubic"
            a = getA();
            inv_d2 = (h.^2 + k.^2 + l.^2) ./ a^2;

        case "tetragonal"
            a = getA(); c = getC();
            inv_d2 = (h.^2 + k.^2) / a^2 + l.^2 / c^2;

        case "orthorhombic"
            a = getA(); b = getB(); c = getC();
            inv_d2 = h.^2 / a^2 + k.^2 / b^2 + l.^2 / c^2;

        case "hexagonal"
            a = getA(); c = getC();
            inv_d2 = (4/3) * (h.^2 + k.^2 + h.*k) / a^2 + l.^2 / c^2;

        case "rhombohedral"
            a = getA();
            % d² = a²·(1 − 3cos²α + 2cos³α) / [(h²+k²+l²)sin²α + 2(hk+kl+lh)(cos²α − cosα)]
            num = a^2 * (1 - 3*cos(alpha)^2 + 2*cos(alpha)^3);
            den = (h.^2 + k.^2 + l.^2) * sin(alpha)^2 + ...
                  2 * (h.*k + k.*l + l.*h) * (cos(alpha)^2 - cos(alpha));
            d2 = num ./ den;
            d = sqrt(d2);
            d(~isfinite(d) | d <= 0) = NaN;
            return

        case "monoclinic"
            a = getA(); b = getB(); c = getC();
            sb = sin(beta); cb = cos(beta);
            inv_d2 = h.^2 ./ (a * sb)^2 + k.^2 / b^2 + l.^2 ./ (c * sb)^2 ...
                    - 2 * h .* l * cb / (a * c * sb^2);

        case "triclinic"
            a = getA(); b = getB(); c = getC();
            ca = cos(alpha); cb = cos(beta); cg = cos(gamma);
            sa = sin(alpha); sb = sin(beta); sg = sin(gamma);

            % Unit cell volume factor:
            %   V² = a²b²c²·(1 − cos²α − cos²β − cos²γ + 2cosα cosβ cosγ)
            volFactor = 1 - ca^2 - cb^2 - cg^2 + 2*ca*cb*cg;

            num = volFactor;
            den = ( (h .* sa / a).^2 + (k .* sb / b).^2 + (l .* sg / c).^2 ) ...
                + ( 2 * h .* k / (a * b) * (ca*cb - cg) ...
                  + 2 * k .* l / (b * c) * (cb*cg - ca) ...
                  + 2 * h .* l / (a * c) * (cg*ca - cb) );
            d2 = num ./ den;
            d = sqrt(d2);
            d(~isfinite(d) | d <= 0) = NaN;
            return

        otherwise
            error('xrdc:lattice:unknownSystem', ...
                'Unknown crystal system: %s', lattice.system);
    end

    d = 1 ./ sqrt(inv_d2);
    d(~isfinite(d) | d <= 0) = NaN;
end

% ---------- helpers ----------

function val = mustHave(lattice, field)
    if ~isfield(lattice, field) || isempty(lattice.(field))
        error('xrdc:lattice:missingField', ...
            'lattice.%s is required for %s system.', field, lattice.system);
    end
    val = lattice.(field);
end

function rad = getAngleRad(lattice, field, defaultDeg)
    if isfield(lattice, field) && ~isempty(lattice.(field))
        rad = deg2rad(lattice.(field));
    else
        rad = deg2rad(defaultDeg);
    end
end

function [h, k, l] = broadcastTriplet(h, k, l)
    % replicate scalars to match any array inputs
    sizes = {size(h), size(k), size(l)};
    nonScalar = cellfun(@(s) prod(s) > 1, sizes);
    if ~any(nonScalar)
        return
    end
    target = sizes{find(nonScalar, 1)};
    if isscalar(h), h = repmat(h, target); end
    if isscalar(k), k = repmat(k, target); end
    if isscalar(l), l = repmat(l, target); end
    assert(isequal(size(h), size(k), size(l)), ...
        'xrdc:lattice:sizeMismatch', 'h, k, l must be broadcastable.');
end
