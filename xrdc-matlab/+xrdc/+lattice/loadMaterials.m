function M = loadMaterials(name)
%LOADMATERIALS  Load the material database (+xrdc/+data/materials.json).
%   M = xrdc.lattice.loadMaterials()        — struct array, all entries
%   e = xrdc.lattice.loadMaterials(name)    — single entry by name/alias
%
%   Lookup is case-insensitive over .name and .aliases. Errors with
%   xrdc:lattice:unknownMaterial (message lists valid names) on no match.
%
%   See also xrdc.lattice.identifyMaterial.

    arguments
        name (1,1) string = ""
    end

    persistent cache
    if isempty(cache)
        jsonPath = fullfile(fileparts(mfilename('fullpath')), ...
            '..', '+data', 'materials.json');
        raw = jsondecode(fileread(jsonPath));
        mats = raw.materials;
        if iscell(mats)   % defensive: heterogeneous fields decode as cell
            mats = [mats{:}];
        end
        cache = mats(:).';
    end

    if name == ""
        M = cache;
        return
    end

    target = lower(name);
    for e = cache
        candidates = lower([string(e.name); string(e.aliases(:))]);
        if any(candidates == target)
            M = e;
            return
        end
    end
    error('xrdc:lattice:unknownMaterial', ...
        'Unknown material "%s". Valid names: %s.', name, ...
        strjoin(string({cache.name}), ', '));
end
