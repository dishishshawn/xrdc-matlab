function f = elasticFactor(material)
%ELASTICFACTOR  Biaxial strain factor f from a material's elastic block.
%   f = xrdc.lattice.elasticFactor(material)
%
%   f relates out-of-plane to in-plane strain as eps_perp = -f * eps_par for
%   a biaxially strained (001) film. Two supported elastic models:
%     elastic.nu          -> f = 2*nu/(1-nu)     (isotropic Poisson)
%     elastic.{c13,c33}   -> f = 2*c13/c33       (tetragonal single-crystal)
%
%   Input
%     material : struct with an .elastic field (e.g. from loadMaterials).
%
%   Errors: xrdc:lattice:noElastic if neither model is present.

    arguments
        material (1,1) struct
    end
    if ~isfield(material, 'elastic') || ~isstruct(material.elastic)
        error('xrdc:lattice:noElastic', ...
            'Material has no elastic block; cannot compute strain factor.');
    end
    e = material.elastic;
    if isfield(e, 'nu') && ~isempty(e.nu)
        f = 2 * e.nu / (1 - e.nu);
    elseif isfield(e, 'c13') && isfield(e, 'c33') && ~isempty(e.c13) && ~isempty(e.c33)
        f = 2 * e.c13 / e.c33;
    else
        error('xrdc:lattice:noElastic', ...
            'elastic block lacks both nu and {c13,c33}.');
    end
end
