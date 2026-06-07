function fig = showSplash()
%SHOWSPLASH  Display the animated XRDC launch splash in a small uifigure.
%
%   fig = xrdc.ui.showSplash() pops a borderless-as-possible launch window
%   showing the self-contained animated splash (resources/splash.html) and
%   returns its uifigure handle. The HTML animates inside an embedded
%   Chromium (uihtml) running in its own process, so the splash keeps
%   playing smoothly while the main app builds on MATLAB's single thread.
%
%   The handle is ALWAYS deletable: if the splash asset can't be located
%   (e.g. it was dropped from a build), an empty handle is returned so the
%   caller can blindly `delete(fig)` without special-casing failure.
%
%   Works both from source and inside the compiled standalone — the HTML is
%   read into memory and passed as HTMLSource, so there are no runtime
%   relative-asset path issues, and the file itself is located by a
%   recursive search under ctfroot (deployed) or the repo root (source).

    fig = gobjects(1, 0);   % deletable no-op if anything below bails

    html = readSplashHtml();
    if html == ""
        return;   % asset missing — caller still gets a valid (empty) handle
    end

    sz = 460;   % square, matches the 1:1 splash artwork
    fig = uifigure( ...
        'Name',     'XRDC', ...
        'Position', centerRect(sz, sz), ...
        'Color',    [0.016 0.024 0.055], ...   % #04060E — no white flash before paint
        'Resize',   'off');
    % Match the artwork's dark background so the (unavoidable) title bar and
    % any pre-paint frame read as part of the splash. Theme is R2025a+.
    try
        fig.Theme = 'dark';
    catch
    end

    uihtml(fig, 'HTMLSource', char(html), 'Position', [0 0 sz sz]);
    drawnow;   % force the splash to paint before we return and build the app
end

% ---------------------------------------------------------------------
function html = readSplashHtml()
%READSPLASHHTML  Load splash.html as a string, or "" if it can't be found.
    html = "";
    if isdeployed
        base = ctfroot;
    else
        here = fileparts(mfilename('fullpath'));   % +xrdc/+ui
        base = fileparts(fileparts(here));         % repo root
    end
    hits = dir(fullfile(base, '**', 'splash.html'));
    if isempty(hits)
        return;
    end
    html = string(fileread(fullfile(hits(1).folder, hits(1).name)));
end

% ---------------------------------------------------------------------
function r = centerRect(w, h)
%CENTERRECT  [x y w h] centering a w-by-h window on the primary screen.
    ss = get(groot, 'ScreenSize');   % [1 1 W H]
    x  = max(1, round((ss(3) - w) / 2));
    y  = max(1, round((ss(4) - h) / 2));
    r  = [x y w h];
end
