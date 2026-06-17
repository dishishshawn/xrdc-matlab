function p = configPath()
%CONFIGPATH  Absolute path of the persisted lab-preset settings file.
%   p = xrdc.config.configPath()
%
%   The settings live as JSON under MATLAB's per-user preferences directory
%   (prefdir), so they persist across sessions and survive a fresh clone of
%   the repo without being committed. The folder is created on demand by
%   xrdc.config.save.

    p = string(fullfile(prefdir, 'xrdc', 'settings.json'));
end
