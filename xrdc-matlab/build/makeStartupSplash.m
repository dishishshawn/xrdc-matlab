function makeStartupSplash()
%MAKESTARTUPSPLASH  Render resources/splash_startup.png from splash_static.html.
%
%   The compiled exe shows a STATIC image while the MATLAB Runtime starts up
%   (before any app code can run, so it can't be the animated HTML splash).
%   buildStandalone passes this PNG as ExecutableSplashScreen. We render it
%   from splash_static.html — a still, text-free unit cell (no wordmark or
%   progress bar) — so the startup image is a clean lead-in to the animated
%   splash rather than a stock MATLAB splash.
%
%   Run this whenever resources/splash_static.html changes, then rebuild:
%     >> addpath build; makeStartupSplash; buildStandalone
%
%   Size (540x540) and background match xrdc.ui.showSplash so the startup
%   image and the animated window overlap seamlessly when both are centered.

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    html = string(fileread(fullfile(repoRoot, 'resources', 'splash_static.html')));

    sz = 540;
    f = uifigure('Name', 'XRDC', 'Position', [200 200 sz sz], ...
        'Color', [0.024 0.035 0.071], 'Resize', 'off');   % ~ #060912 splash bg
    cleanup = onCleanup(@() delete(f));
    try, f.Theme = 'dark'; catch, end
    uihtml(f, 'HTMLSource', char(html), 'Position', [0 0 sz sz]);
    drawnow;
    pause(0.6);   % static art: just let the SVG build + paint

    out = fullfile(repoRoot, 'resources', 'splash_startup.png');
    exportapp(f, out);
    fprintf('Wrote startup splash: %s\n', out);
end
