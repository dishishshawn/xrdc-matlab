function cfg = load(path)
%LOAD  Read persisted lab presets, merged over the factory defaults.
%   cfg = xrdc.config.load()
%   cfg = xrdc.config.load(path)
%
%   Returns xrdc.config.defaults() with any saved overrides applied. Unknown
%   or malformed fields in the file are ignored (defaults kept) so an old or
%   hand-edited settings file can never crash the app, and new default fields
%   added in a later version are picked up automatically. A corrupt JSON file
%   degrades to defaults with a warning rather than an error.
%
%   path defaults to xrdc.config.configPath().

    arguments
        path (1,1) string = xrdc.config.configPath()
    end

    cfg = xrdc.config.defaults();
    if ~isfile(path)
        return
    end

    try
        raw    = fileread(char(path));
        loaded = jsondecode(raw);
    catch ME
        warning('xrdc:config:badFile', ...
            'Could not parse settings file %s (%s); using defaults.', ...
            path, ME.message);
        return
    end

    if ~isstruct(loaded)
        warning('xrdc:config:badFile', ...
            'Settings file %s is not a JSON object; using defaults.', path);
        return
    end

    % Merge only known fields, coercing each to the default's type so a
    % stray string/number in the file can't change a field's class.
    for fn = fieldnames(cfg).'
        name = fn{1};
        if ~isfield(loaded, name), continue, end
        cfg.(name) = coerce(cfg.(name), loaded.(name));
    end
end

% =====================================================================
function out = coerce(template, value)
%COERCE  Cast a loaded value to the class of the default it overrides.
    if islogical(template)
        out = logical(value);
    elseif isstring(template) || ischar(template)
        out = string(value);
    elseif isnumeric(template)
        v = double(value);
        if isscalar(v) && isfinite(v)
            out = v;
        else
            out = template;   % reject non-scalar / non-finite numerics
        end
    else
        out = value;
    end
end
