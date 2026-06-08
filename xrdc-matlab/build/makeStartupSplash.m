function makeStartupSplash()
%MAKESTARTUPSPLASH  Render resources/splash_startup.png from splash.html.
%
%   The compiled exe shows a STATIC image while the MATLAB Runtime starts up
%   (before any app code can run, so it can't be the animated HTML splash).
%   buildStandalone passes this PNG as ExecutableSplashScreen so that static
%   startup image is a still frame of the animated splash — the two then read
%   as one continuous launch screen instead of a stock MATLAB splash.
%
%   Run this whenever resources/splash.html changes, then rebuild:
%     >> addpath build; makeStartupSplash; buildStandalone
%
%   Size (460x460) and background match xrdc.ui.showSplash so the startup
%   image and the animated window overlap seamlessly when both are centered.

    repoRoot = fileparts(fileparts(mfilename('fullpath')));
    html = string(fileread(fullfile(repoRoot, 'resources', 'splash.html')));

    sz = 540;
    f = uifigure('Name', 'XRDC', 'Position', [200 200 sz sz], ...
        'Color', [0.024 0.035 0.071], 'Resize', 'off');   % ~ #060912 splash bg
    cleanup = onCleanup(@() delete(f));
    try, f.Theme = 'dark'; catch, end
    uihtml(f, 'HTMLSource', char(html), 'Position', [0 0 sz sz]);
    drawnow;
    pause(1.3);   % let the launch + rise finish and the cell settle into spin

    out = fullfile(repoRoot, 'resources', 'splash_startup.png');
    exportapp(f, out);
    fprintf('Wrote startup splash: %s\n', out);
end
