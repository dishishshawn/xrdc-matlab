function save(cfg, path)
%SAVE  Persist lab presets to JSON.
%   xrdc.config.save(cfg)
%   xrdc.config.save(cfg, path)
%
%   Writes cfg (a struct in the xrdc.config.defaults() shape) to the user's
%   settings file, creating the containing folder on demand. Only known
%   fields are written, so an over-stuffed struct can't pollute the file;
%   the values are pretty-printed for hand-editing.
%
%   path defaults to xrdc.config.configPath().

    arguments
        cfg  (1,1) struct
        path (1,1) string = xrdc.config.configPath()
    end

    % Keep only the canonical fields, in a stable order, so the file stays
    % a clean record of the documented settings.
    def = xrdc.config.defaults();
    out = struct();
    for fn = fieldnames(def).'
        name = fn{1};
        if isfield(cfg, name)
            out.(name) = cfg.(name);
        else
            out.(name) = def.(name);
        end
    end

    folder = fileparts(path);
    if strlength(folder) > 0 && ~isfolder(folder)
        [ok, msg] = mkdir(char(folder));
        if ~ok
            error('xrdc:config:writeFailed', ...
                'Could not create settings folder %s: %s', folder, msg);
        end
    end

    try
        json = jsonencode(out, 'PrettyPrint', true);
    catch
        json = jsonencode(out);   % older releases lack PrettyPrint
    end

    fid = fopen(char(path), 'w');
    if fid < 0
        error('xrdc:config:writeFailed', ...
            'Could not open settings file for writing: %s', path);
    end
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, json, 'char');
end
