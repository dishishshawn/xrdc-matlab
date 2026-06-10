function xrdcApp()
%XRDCAPP  GUI for the xrdc-matlab toolkit — load, analyze, export XRD scans.
%
%   Usage:
%     >> xrdcApp
%
%   Opens an interactive window. Click "Load Scan..." to pick a file;
%   the app detects scan type (rocking curve / θ-2θ / phi / XRR / RSM)
%   from the file contents and filename, runs the appropriate analysis,
%   and shows a live plot. Tweak the parameters on the left, click
%   "Export" to save a 600 dpi publication-quality image.
%
%   Designed so that a lab member with zero MATLAB knowledge can process
%   a new Rigaku / PANalytical file and get a paper-ready figure without
%   editing any code.

    thisDir = fileparts(mfilename('fullpath'));
    addpath(thisDir);

    % Animated launch splash. Shown first; the main window is built hidden
    % and revealed at the end so the splash stays on top while loading.
    splashFig = xrdc.ui.showSplash();
    splashStart = tic;

    T = appTheme();

    fig = uifigure('Name', 'XRDC Scan Analyzer', 'Position', [100 100 1240 780], ...
        'Visible', 'off', 'Color', T.bg);
    % Keep a light THEME so the plot axes default to white (publication-ready);
    % the dark crystalline "chrome" below is layered on with explicit colours,
    % so the OS dark mode can't wash out the plot. Theme is R2025a+; guard so
    % older releases keep their default.
    try
        fig.Theme = 'light';
    catch
    end

    grid = uigridlayout(fig, [4 2]);
    grid.RowHeight     = {54, 40, 26, '1x'};
    grid.ColumnWidth   = {310, '1x'};
    grid.RowSpacing    = 8;
    grid.ColumnSpacing = 8;
    grid.Padding       = [10 10 10 10];
    grid.BackgroundColor = T.bg;

    % Branded header: crystal wordmark in the splash palette
    header = uigridlayout(grid, [1 1]);
    header.Layout.Row = 1; header.Layout.Column = [1 2];
    header.Padding = [4 2 4 2];
    header.BackgroundColor = T.bg;
    uilabel(header, 'Text', '◆ XRDC', ...
        'FontName', T.font, 'FontSize', 24, 'FontWeight', 'bold', ...
        'FontColor', T.gold, 'VerticalAlignment', 'center');

    % Top bar: load (primary, gold) + export + customize (secondary)
    topBar = uigridlayout(grid, [1 4]);
    topBar.Layout.Row = 2; topBar.Layout.Column = [1 2];
    topBar.ColumnWidth = {150, 160, 170, '1x'};
    topBar.ColumnSpacing = 8; topBar.Padding = [0 0 0 0];
    topBar.BackgroundColor = T.bg;
    uibutton(topBar, 'Text', 'Load Scan...', ...
        'FontName', T.font, 'FontSize', 13, 'FontWeight', 'bold', ...
        'BackgroundColor', T.gold, 'FontColor', T.ink, ...
        'ButtonPushedFcn', @(~,~) onLoadScan(fig));
    exportBtn = uibutton(topBar, 'Text', 'Export 600 dpi...', ...
        'FontName', T.font, 'FontSize', 13, 'Enable', 'off', ...
        'BackgroundColor', T.btn, 'FontColor', T.text, ...
        'ButtonPushedFcn', @(~,~) onExport(fig));
    customizeBtn = uibutton(topBar, 'Text', 'Customize Plot...', ...
        'FontName', T.font, 'FontSize', 13, 'Enable', 'off', ...
        'BackgroundColor', T.btn, 'FontColor', T.text, ...
        'ButtonPushedFcn', @(~,~) onCustomizePlot(fig));

    % Info strip
    infoLbl = uilabel(grid, 'Text', '  No scan loaded. Click "Load Scan..." to begin.', ...
        'FontName', T.font, 'FontSize', 12, 'FontColor', T.textDim, ...
        'HorizontalAlignment', 'left');
    infoLbl.Layout.Row = 3; infoLbl.Layout.Column = [1 2];

    % Left: analysis panel
    leftPanel = uipanel(grid, 'Title', 'Analysis', ...
        'FontName', T.font, 'FontSize', 13, 'FontWeight', 'bold', ...
        'BackgroundColor', T.panel, 'ForegroundColor', T.gold, ...
        'BorderColor', T.edge);
    leftPanel.Layout.Row = 4; leftPanel.Layout.Column = 1;

    % Right: plot preview (uiaxes stays white → publication figure on a card)
    plotPanel = uipanel(grid, 'Title', 'Preview', ...
        'FontName', T.font, 'FontSize', 13, 'FontWeight', 'bold', ...
        'BackgroundColor', T.panel, 'ForegroundColor', T.gold, ...
        'BorderColor', T.edge);
    plotPanel.Layout.Row = 4; plotPanel.Layout.Column = 2;
    plotGrid = uigridlayout(plotPanel, [1 1]);
    % WHITE plot card: the axes box is white, but the title/axis-label margin
    % shows this background — keep it white so black title/labels stay visible
    % (and the on-screen preview matches the exported figure). The dark panel
    % frame around it makes it read as a card.
    plotGrid.Padding = [6 6 6 6]; plotGrid.BackgroundColor = [1 1 1];
    ax = uiaxes(plotGrid);

    % Store state on the figure so callbacks can share it
    st = struct();
    st.scan        = [];
    st.rsmScans    = [];
    st.overlayScan = [];      % 2nd RC for the film-vs-substrate overlay
    st.overlayPath = '';
    st.filePath    = '';
    st.detectedType = "";
    st.params      = struct();
    st.style       = defaultStyle();
    st.ax          = ax;
    st.infoLbl     = infoLbl;
    st.exportBtn   = exportBtn;
    st.customizeBtn = customizeBtn;
    st.leftPanel   = leftPanel;
    st.resultsArea = [];
    fig.UserData = st;

    placeholder(ax);

    % Reveal the finished app, then dismiss the splash — but keep the splash
    % up for a minimum so its launch animation is actually seen on fast loads.
    minSplashSecs = 1.8;
    pause(max(0, minSplashSecs - toc(splashStart)));
    fig.Visible = 'on';
    drawnow;
    delete(splashFig);
end

% =====================================================================
% Callbacks
% =====================================================================
function onLoadScan(fig)
    [file, path] = uigetfile({ ...
        '*.txt;*.xrdml;*.xrdc;*.x00', 'XRD scan files (*.txt, *.xrdml)'; ...
        '*.*', 'All files (*.*)'}, ...
        'Select an XRD scan');
    if isequal(file, 0), return, end
    fullPath = fullfile(path, file);

    st = fig.UserData;
    st.filePath = fullPath;
    st.params   = struct();   % reset per-scan parameters
    st.style    = defaultStyle();   % reset plot-style overrides per scan
    st.rsmScans = [];
    st.overlayScan = [];      % drop any overlay from a previous scan
    st.overlayPath = '';

    dlg = uiprogressdlg(fig, 'Title', 'Loading scan', ...
        'Message', sprintf('Reading %s ...', file), ...
        'Indeterminate', 'on', 'Cancelable', 'off');
    closeDlg = onCleanup(@() close(dlg));

    try
        fnLower = lower(string(file));
        if endsWith(fnLower, '.xrdml') && contains(fnLower, 'rsm')
            st.rsmScans    = xrdc.rsm.loadAreaScan({fullPath});
            st.scan        = st.rsmScans(1);   % representative
            st.detectedType = "rsm";
        else
            st.scan        = xrdc.io.readScan(fullPath);
            st.detectedType = detectScanType(st.scan, fnLower);
        end
    catch ME
        uialert(fig, sprintf('Failed to load the file:\n\n%s', ME.message), ...
            'Load error', 'Icon', 'error');
        return
    end

    st.infoLbl.Text = sprintf('  %s     |     detected: %s     |     %d points', ...
        file, upper(char(st.detectedType)), numel(st.scan.twoTheta));
    st.infoLbl.FontColor = appTheme().text;   % icy text on the dark chrome
    st.exportBtn.Enable = 'on';
    st.customizeBtn.Enable = 'on';
    fig.UserData = st;

    dlg.Message = 'Running analysis ...';
    drawnow
    buildAnalysisPanel(fig);
    runAnalysis(fig);
end

function onExport(fig)
    st = fig.UserData;
    if isempty(st.scan), return, end

    [~, stem] = fileparts(st.filePath);
    prefix = prefixForType(st.detectedType);
    defaultName = sprintf('%s_%s.png', prefix, stem);

    [file, path] = uiputfile({ ...
        '*.png', 'PNG image 600 dpi (*.png)'; ...
        '*.pdf', 'PDF vector (*.pdf)'; ...
        '*.svg', 'SVG vector (*.svg)'}, ...
        'Export figure', defaultName);
    if isequal(file, 0), return, end

    target = fullfile(path, file);
    W = styleNum(st.style, 'exportW');
    H = styleNum(st.style, 'exportH');
    try
        if ~isnan(W) && ~isnan(H) && W > 0 && H > 0
            exportSized(st.ax, target, W, H);   % honor Customize Plot size
        else
            exportgraphics(st.ax, target, 'Resolution', 600);
        end
        uialert(fig, sprintf('Saved:\n%s', target), 'Export', 'Icon', 'success');
    catch ME
        uialert(fig, sprintf('Export failed:\n%s', ME.message), 'Error', 'Icon', 'error');
    end
end

function exportSized(ax, target, widthIn, heightIn)
%EXPORTSIZED  Export an axes at a specific physical size (inches @ 600 dpi).
%   copyobj the axes into a temporary uifigure sized to the requested
%   inches so the exported figure matches a journal column width.
    tmp = uifigure('Visible', 'off', ...
        'Position', [0 0 round(widthIn * 96) round(heightIn * 96)]);
    cleanup = onCleanup(@() delete(tmp));
    tg = uigridlayout(tmp, [1 1], 'Padding', [2 2 2 2]);
    ax2 = copyobj(ax, tg); %#ok<NASGU>
    exportgraphics(ax2, target, 'Resolution', 600);
end

function onParamChange(fig, name, value)
    st = fig.UserData;
    st.params.(name) = value;
    fig.UserData = st;

    dlg = uiprogressdlg(fig, 'Title', 'Updating', ...
        'Message', 'Re-running analysis ...', ...
        'Indeterminate', 'on', 'Cancelable', 'off');
    closeDlg = onCleanup(@() close(dlg));
    drawnow

    runAnalysis(fig);
end

function onRcOverlayToggle(fig, src)
%ONRCOVERLAYTOGGLE  Checkbox handler for the film-vs-substrate RC overlay.
%   Turning it on prompts for the complementary RC (asks for the film when a
%   substrate is loaded, and vice-versa, inferred from the loaded filename).
%   Cancelling the file dialog reverts the checkbox. Turning it off drops the
%   overlay. Either way the rocking-curve analysis re-renders.
    st = fig.UserData;

    if src.Value
        [~, baseName] = fileparts(st.filePath);
        want = complementaryRcRole(baseName);   % "film" / "substrate" / "corresponding"
        [file, path] = uigetfile({ ...
            '*.txt;*.xrdml;*.xrdc;*.x00', 'XRD scan files (*.txt, *.xrdml)'; ...
            '*.*', 'All files (*.*)'}, ...
            sprintf('Select the %s rocking curve to overlay', want));
        if isequal(file, 0)
            src.Value = false;          % user cancelled → leave overlay off
            return
        end
        try
            ov = xrdc.io.readScan(fullfile(path, file));
        catch ME
            uialert(fig, sprintf('Could not read the overlay RC:\n\n%s', ME.message), ...
                'Overlay error', 'Icon', 'error');
            src.Value = false;
            return
        end
        st.overlayScan      = ov;
        st.overlayPath      = fullfile(path, file);
        st.params.rcOverlay = true;
    else
        st.overlayScan      = [];
        st.overlayPath      = '';
        st.params.rcOverlay = false;
    end
    fig.UserData = st;

    dlg = uiprogressdlg(fig, 'Title', 'Updating', ...
        'Message', 'Rendering overlay ...', 'Indeterminate', 'on', 'Cancelable', 'off');
    closeDlg = onCleanup(@() close(dlg));
    drawnow
    runAnalysis(fig);
end

% =====================================================================
% Panel construction per scan type
% =====================================================================
function buildAnalysisPanel(fig)
    st = fig.UserData;
    T  = appTheme();
    delete(st.leftPanel.Children);

    g = uigridlayout(st.leftPanel, [14 2]);
    g.RowHeight   = [repmat({26}, 1, 6), {'1x'}, repmat({26}, 1, 7)];
    g.ColumnWidth = {110, '1x'};
    g.RowSpacing  = 4;
    g.Padding     = [6 6 6 6];
    g.BackgroundColor = T.panel;

    row = 1;
    t = char(lower(string(st.detectedType)));
    switch t
        case 'omega'
            row = addEdit (g, row, 'Fit window (°)', '0.5', @(v) onParamChange(fig, 'fitWindow', v));
            row = addDrop (g, row, 'Shape', {'gauss','lorentz','pseudoVoigt'}, 'gauss', ...
                                                      @(v) onParamChange(fig, 'shape',     v));
            showFit = isfield(st.params, 'showFit') && st.params.showFit;
            row = addCheck(g, row, 'Show fit curve', showFit, ...
                                                      @(src) onParamChange(fig, 'showFit', src.Value));
            overlayOn = isfield(st.params, 'rcOverlay') && st.params.rcOverlay;
            row = addCheck(g, row, 'Overlay 2nd RC (film/sub)', overlayOn, ...
                                                      @(src) onRcOverlayToggle(fig, src));
        case 'twothetaomega'
            row = addEdit (g, row, 'Min prom (%)',   '5',   @(v) onParamChange(fig, 'promPct',   v));
            subs = identifiableSubstrates();
            row = addDrop (g, row, 'Substrate', subs, subs{1}, ...
                                                      @(v) onSubstrateChange(fig, v));
            row = addIdentifyButton(g, row, fig); %#ok<NASGU>
        case 'xrr'
            row = addEdit (g, row, 'Fringe 2θ min',  '0',   @(v) onParamChange(fig, 'xrrMin',    v));
            row = addEdit (g, row, 'Fringe 2θ max',  '5.0', @(v) onParamChange(fig, 'xrrMax',    v));
            row = addEdit (g, row, 'Min prom (%)',   '1.5', @(v) onParamChange(fig, 'xrrProm',   v));
        case 'phi'
            row = addEdit (g, row, 'Noise σ (k)',    '6',   @(v) onParamChange(fig, 'noiseSigmas', v));
        case 'rsm'
            row = addDrop (g, row, 'Colormap', {'turbo','parula','jet'}, 'turbo', ...
                                                      @(v) onParamChange(fig, 'colormap', v));
            row = addEdit (g, row, 'Imin (counts)',  '10',  @(v) onParamChange(fig, 'imin',      v));
            row = addEdit (g, row, 'Imax (counts)',  '1e5', @(v) onParamChange(fig, 'imax',      v));
            row = addEdit (g, row, 'Contours',       '40',  @(v) onParamChange(fig, 'nContours', v));
    end

    % Results area fills the rest
    hdr = uilabel(g, 'Text', 'Results', 'FontWeight', 'bold', ...
        'FontName', T.font, 'FontColor', T.gold);
    hdr.Layout.Row = 8; hdr.Layout.Column = [1 2];
    ta = uitextarea(g, 'Editable', 'off', 'FontName', T.mono, 'FontSize', 11, ...
        'BackgroundColor', T.panel2, 'FontColor', T.text);
    ta.Layout.Row = [9 14]; ta.Layout.Column = [1 2];
    st.resultsArea = ta;
    fig.UserData = st;
end

function row = addEdit(g, row, label, default, cb)
    T = appTheme();
    lbl = uilabel(g, 'Text', label, 'FontName', T.font, 'FontColor', T.text);
    lbl.Layout.Row = row; lbl.Layout.Column = 1;
    ef = uieditfield(g, 'Value', default, ...
        'ValueChangedFcn', @(src, ~) cb(src.Value));
    ef.Layout.Row = row; ef.Layout.Column = 2;
    row = row + 1;
end

function row = addDrop(g, row, label, items, default, cb)
    T = appTheme();
    lbl = uilabel(g, 'Text', label, 'FontName', T.font, 'FontColor', T.text);
    lbl.Layout.Row = row; lbl.Layout.Column = 1;
    dd = uidropdown(g, 'Items', items, 'Value', default, ...
        'ValueChangedFcn', @(src, ~) cb(src.Value));
    dd.Layout.Row = row; dd.Layout.Column = 2;
    row = row + 1;
end

function row = addCheck(g, row, label, value, cb)
    T = appTheme();
    chk = uicheckbox(g, 'Text', label, 'Value', value, ...
        'FontName', T.font, 'FontColor', T.text, ...
        'ValueChangedFcn', @(src, ~) cb(src));
    chk.Layout.Row = row; chk.Layout.Column = [1 2];
    row = row + 1;
end

function subs = identifiableSubstrates()
%IDENTIFIABLESUBSTRATES  Dropdown items: materials usable as a substrate.
    M = xrdc.lattice.loadMaterials();
    ok = ismember(string({M.role}), ["substrate", "both"]);
    subs = cellstr(string({M(ok).name}));
end

function row = addIdentifyButton(g, row, fig)
    T = appTheme();
    b = uibutton(g, 'Text', 'Identify Material', ...
        'FontName', T.font, 'FontSize', 12, ...
        'BackgroundColor', T.btn, 'FontColor', T.text, ...
        'ButtonPushedFcn', @(~,~) onIdentifyMaterial(fig));
    b.Layout.Row = row; b.Layout.Column = [1 2];
    row = row + 1;
end

function T = appTheme()
%APPTHEME  Dark crystalline UI palette, derived from the launch splash.
%   The plot axes deliberately stay white (publication-ready); only the
%   surrounding "chrome" uses these colours. Chakra Petch (the splash font)
%   is a web font unavailable to native uicomponents, so the chrome uses the
%   system 'Segoe UI'; the splash itself keeps its web fonts via uihtml.
    % Softened palette: muted slate-navy (not near-black) and a warm,
    % de-neoned gold, for an easier-on-the-eyes look.
    T.bg      = [0.106 0.129 0.184];   % #1B212F  window background
    T.panel   = [0.141 0.169 0.231];   % #242B3B  panels
    T.panel2  = [0.122 0.149 0.204];   % #1F2634  deeper, results readout
    T.btn     = [0.196 0.235 0.318];   % #323C51  secondary buttons
    T.edge    = [0.247 0.290 0.380];   % #3F4A61  soft panel borders
    T.text    = [0.851 0.886 0.933];   % #D9E2EE  soft icy text
    T.textDim = [0.557 0.604 0.702];   % #8E9AB3  dimmed labels
    T.gold    = [0.882 0.741 0.420];   % #E1BD6B  soft warm gold accent
    T.blue    = [0.490 0.651 0.851];   % #7DA6D9  soft blue accent
    T.ink     = [0.106 0.129 0.184];   % dark text on the gold button
    T.font    = 'Segoe UI';
    T.mono    = 'Consolas';
end

% =====================================================================
% Analysis dispatchers
% =====================================================================
function runAnalysis(fig)
    st = fig.UserData;
    cla(st.ax); reset(st.ax);
    % Dispatch on a char vector to avoid any string/char switch quirks.
    t = char(lower(string(st.detectedType)));
    try
        switch t
            case 'omega',         runRockingCurve(fig);
            case 'twothetaomega', runThetaTwoTheta(fig);
            case 'xrr',           runXRR(fig);
            case 'phi',           runPhiScan(fig);
            case 'rsm',           runRSM(fig);
            otherwise,            plotBasic(fig);
        end
    catch ME
        uialert(fig, sprintf('Analysis error:\n\n%s', ME.message), ...
            'Analysis error', 'Icon', 'error');
    end

    % Apply the user's optional plot-style overrides on top of the
    % per-analysis defaults (no-op when nothing is overridden).
    try
        s = fig.UserData;
        xrdc.plot.applyStyle(s.ax, s.style);
    catch
    end
end

function t = detectScanType(scan, fnLower)
%DETECTSCANTYPE  Pick the right analysis route from filename + scan.scanType.
%
%   Filename keywords win over scanType because Rigaku labels everything
%   "twoThetaOmega" even when the scan is really XRR or an RC, and the
%   filename is what the user actually recognises.

    if contains(fnLower, "xrr")
        t = "xrr";                   return
    end
    if contains(fnLower, "rsm")
        t = "rsm";                   return
    end
    if contains(fnLower, ["phi", "φ"])
        t = "phi";                   return
    end
    if contains(fnLower, [" rc", "_rc", "rocking"])
        t = "omega";                 return
    end
    if contains(fnLower, ["2theta", "th2th", "2th", "theta-omega"])
        t = "twothetaomega";         return
    end
    % Fall back to the loader's own inference
    t = lower(string(scan.scanType));
    if t == "twothetaomega"   % leave as-is
        return
    elseif t == "omega" || t == "phi" || t == "area"
        return
    else
        t = "twothetaomega";         % safe default: treat as θ-2θ
    end
end

function runRockingCurve(fig)
    st   = fig.UserData;
    scan = st.scan;    ax = st.ax;

    w     = getNum(st.params, 'fitWindow', 0.5);
    shape = getStr(st.params, 'shape',    'gauss');

    % Overlay mode: compare a film RC against its substrate RC on one set of
    % normalized, peak-centered axes — the standard deposition-quality check.
    overlayOn = isfield(st, 'overlayScan') && ~isempty(st.overlayScan) ...
        && isfield(st.params, 'rcOverlay') && st.params.rcOverlay;
    if overlayOn
        runRockingCurveOverlay(fig, w, shape);
        return
    end

    pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts) * 0.05);
    semilogy(ax, scan.twoTheta, max(scan.counts, 1), '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5);
    stylePubAxes(ax, '\omega (°)', 'Counts', '');

    if isempty(pk)
        writeResults(fig, {'No peak detected — try lowering threshold.'});
        return
    end

    [~, idx] = max([pk.counts]); pkMain = pk(idx);
    window = [pkMain.twoTheta - w, pkMain.twoTheta + w];

    fit = xrdc.peaks.fitPeak(scan, window, 'Shape', string(shape));

    % Cross-check FWHM under the other profile shapes. The shape choice can
    % shift FWHM by ~15-20% (e.g. a Gaussian fit reads wider than a Lorentzian
    % on the same curve), so surface it rather than let it hide behind the
    % default — lets the user reconcile against externally-reported values.
    xcheck = strings(0, 1);
    for s = ["lorentz", "gauss", "pseudoVoigt"]
        if s == string(shape), continue, end
        try
            f2 = xrdc.peaks.fitPeak(scan, window, 'Shape', s);
            xcheck(end+1, 1) = sprintf('%s %.1f″', s, f2.fwhm * 3600); %#ok<AGROW>
        catch
            xcheck(end+1, 1) = sprintf('%s n/a', s); %#ok<AGROW>
        end
    end

    % The fitted profile is intentionally NOT drawn — only the raw RC and its
    % peak marker. The FWHM (from the fit) is still reported in the title/results.
    hold(ax, 'on');
    plot(ax, pkMain.twoTheta, pkMain.counts, 'o', ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', [1 0.8 0.2], 'MarkerSize', 9);
    hold(ax, 'off');
    title(ax, sprintf('Rocking curve — FWHM = %.1f arcsec', fit.fwhm * 3600));

    writeResults(fig, { ...
        sprintf('Peak ω₀     = %.4f °', fit.twoTheta), ...
        sprintf('FWHM        = %.4f ° (%.1f arcsec)', fit.fwhm, fit.fwhm * 3600), ...
        sprintf('Amplitude   = %.2g counts', fit.amplitude), ...
        sprintf('R²          = %.4f', fit.rSquared), ...
        sprintf('Shape       = %s', fit.shape), ...
        sprintf('Other shapes: %s', strjoin(cellstr(xcheck), '   ')), ...
        '', ...
        sprintf('Fit window: ±%.3f° around peak', w)});
end

function runRockingCurveOverlay(fig, w, shape)
%RUNROCKINGCURVEOVERLAY  Film vs. substrate rocking curves on shared axes.
%   Each curve is normalized to its own peak and centered on its own fitted
%   ω₀, so the comparison is purely of mosaic width (FWHM). Substrate renders
%   black, film red (matching the lab's published figure convention); the
%   legend carries both FWHMs and the results panel the FWHM ratio.
    st = fig.UserData;  ax = st.ax;

    [aFit, aOk] = rcFit(st.scan,        w, shape);
    [bFit, bOk] = rcFit(st.overlayScan, w, shape);
    if ~aOk || ~bOk
        writeResults(fig, {'Could not fit a peak in one of the rocking curves.', ...
            'Try widening the fit window or check the overlay file.'});
        return
    end

    [~, aName] = fileparts(st.filePath);
    [~, bName] = fileparts(st.overlayPath);
    [aLabel, aColor] = rcRoleStyle(aName, 1);
    [bLabel, bColor] = rcRoleStyle(bName, 2);

    cla(ax); hold(ax, 'on');
    plotRcNormalized(ax, st.scan,        aFit.twoTheta, aColor);
    plotRcNormalized(ax, st.overlayScan, bFit.twoTheta, bColor);
    hold(ax, 'off');
    set(ax, 'YScale', 'log');
    stylePubAxes(ax, '\Delta\omega (°)', 'Normalized intensity', 'Rocking-curve overlay');

    lg = legend(ax, { ...
        sprintf('%s — FWHM %.1f″', aLabel, aFit.fwhm * 3600), ...
        sprintf('%s — FWHM %.1f″', bLabel, bFit.fwhm * 3600)}, ...
        'Location', 'northeast');
    lg.TextColor = [0 0 0]; lg.Box = 'on';

    ratio = max(aFit.fwhm, bFit.fwhm) / max(min(aFit.fwhm, bFit.fwhm), eps);
    writeResults(fig, { ...
        'RC overlay (each normalized & centered on its peak):', '', ...
        sprintf('%-10s  ω₀ = %8.4f°   FWHM = %.4f° (%.1f″)', ...
            aLabel, aFit.twoTheta, aFit.fwhm, aFit.fwhm * 3600), ...
        sprintf('%-10s  ω₀ = %8.4f°   FWHM = %.4f° (%.1f″)', ...
            bLabel, bFit.twoTheta, bFit.fwhm, bFit.fwhm * 3600), ...
        '', ...
        sprintf('FWHM ratio (broad/narrow) = %.2f', ratio), ...
        sprintf('Peak offset Δω₀          = %.4f°', aFit.twoTheta - bFit.twoTheta), ...
        '', ...
        sprintf('Base    : %s', aName), ...
        sprintf('Overlay : %s', bName), ...
        sprintf('Shape: %s   Fit window: ±%.3f°', shape, w)});
end

function [fit, ok] = rcFit(scan, w, shape)
%RCFIT  Fit the dominant rocking-curve peak; ok=false if none is found.
    ok  = false;
    fit = struct('twoTheta', NaN, 'fwhm', NaN);
    pk  = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts) * 0.05);
    if isempty(pk), return, end
    [~, idx] = max([pk.counts]);
    window   = [pk(idx).twoTheta - w, pk(idx).twoTheta + w];
    try
        fit = xrdc.peaks.fitPeak(scan, window, 'Shape', string(shape));
        ok  = true;
    catch
        ok  = false;
    end
end

function plotRcNormalized(ax, scan, omega0, color)
%PLOTRCNORMALIZED  One RC, peak-normalized and centered on omega0 (log-safe).
    dw = scan.twoTheta - omega0;
    yn = double(scan.counts) / max(double(scan.counts));
    yn(yn <= 0) = NaN;                       % keep the log axis clean
    plot(ax, dw, yn, '-', 'Color', color, 'LineWidth', 1.4);
end

function want = complementaryRcRole(name)
%COMPLEMENTARYRCROLE  What to ask the user to import, given the loaded file.
    n = lower(string(name));
    if contains(n, "film")
        want = 'substrate';
    elseif contains(n, "sub")
        want = 'film';
    else
        want = 'corresponding';
    end
end

function [label, color] = rcRoleStyle(name, ordinal)
%RCROLESTYLE  Map a filename to a (label, colour): substrate→black, film→red.
%   Falls back to a generic "RC <ordinal>" with a distinct colour when the
%   role can't be inferred from the filename.
    n = lower(string(name));
    if contains(n, "film")
        label = "Film";       color = [0.85 0.20 0.20];
    elseif contains(n, "sub")
        label = "Substrate";  color = [0.00 0.00 0.00];
    else
        label = sprintf("RC %d", ordinal);
        if ordinal == 1, color = [0.10 0.40 0.80]; else, color = [0.85 0.20 0.20]; end
    end
end

function runThetaTwoTheta(fig)
    st   = fig.UserData;
    scan = st.scan; ax = st.ax;

    promPct = getNum(st.params, 'promPct', 5);
    pk = xrdc.peaks.findPeaks(scan, ...
        'MinProminence', max(scan.counts) * promPct / 100);

    semilogy(ax, scan.twoTheta, max(scan.counts, 1), '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5);
    hold(ax, 'on');
    if ~isempty(pk)
        plot(ax, [pk.twoTheta], [pk.counts], 'v', ...
            'MarkerFaceColor', [0.85 0.2 0.2], 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    hold(ax, 'off');
    stylePubAxes(ax, '2\theta (°)', 'Counts', sprintf('θ-2θ scan — %d peaks', numel(pk)));

    lines = {sprintf('Detected %d peak(s):', numel(pk)), ''};
    for i = 1:numel(pk)
        lines{end+1} = sprintf('  %2d.  2θ = %7.3f°    I = %7.0f    FWHM = %.3f°', ...
            i, pk(i).twoTheta, pk(i).counts, pk(i).fwhm); %#ok<AGROW>
    end
    writeResults(fig, lines);
end

function onIdentifyMaterial(fig)
%ONIDENTIFYMATERIAL  Run material ID on the current theta-2theta peaks.
    st = fig.UserData;
    if isempty(st.scan), return, end
    promPct = getNum(st.params, 'promPct', 5);
    pk = xrdc.peaks.findPeaks(st.scan, ...
        'MinProminence', max(st.scan.counts) * promPct / 100);
    if isempty(pk)
        uialert(fig, 'No peaks detected - lower "Min prom (%)" and retry.', ...
            'Identify material', 'Icon', 'warning');
        return
    end
    sub = getStr(st.params, 'substrate', '');
    if isempty(sub)
        s = identifiableSubstrates();
        sub = s{1};
    end
    lambda = 1.5406;
    if isfield(st.scan, 'lambda') && ~isempty(st.scan.lambda) ...
            && isfinite(st.scan.lambda)
        lambda = st.scan.lambda;
    end
    try
        R = xrdc.lattice.identifyMaterial(pk, lambda, Substrate=string(sub));
    catch ME
        uialert(fig, sprintf('Identification failed:\n\n%s', ME.message), ...
            'Identify material', 'Icon', 'error');
        return
    end
    annotateIdentification(st.ax, R);
    writeResults(fig, identificationReport(R));
    st = fig.UserData;
    st.params.identified = true;
    fig.UserData = st;
end

function onSubstrateChange(fig, v)
%ONSUBSTRATECHANGE  Substrate dropdown handler.
%   Stores the new substrate; if an identification overlay is currently
%   displayed, re-renders the base plot and re-runs identification so the
%   annotations reflect the new substrate. Without an identification the
%   base plot does not depend on the substrate, so nothing re-renders.
    st = fig.UserData;
    st.params.substrate = v;
    fig.UserData = st;
    if isfield(st.params, 'identified') && st.params.identified
        runAnalysis(fig);
        onIdentifyMaterial(fig);
    end
end

function annotateIdentification(ax, R)
%ANNOTATEIDENTIFICATION  Material + (00l) labels above identified peaks.
    hold(ax, 'on');
    yl = ylim(ax);
    for i = 1:numel(R.substrate.twoTheta)
        labelPeak(ax, R.substrate.twoTheta(i), yl, ...
            sprintf('%s', shortName(R.substrate.name)), [0 0 0]);
    end
    for i = 1:height(R.series)
        name = R.series.bestMatch(i);
        if name == "", name = "?"; end
        tts = seriesRow(R.series.twoTheta, i);  ords = seriesRow(R.series.orders, i);
        for j = 1:numel(tts)
            labelPeak(ax, tts(j), yl, ...
                sprintf('%s (00%d)', shortName(name), ords(j)), ...
                [0.65 0.15 0.15]);
        end
    end
    hold(ax, 'off');
end

function labelPeak(ax, tt, yl, txt, color)
    text(ax, tt, yl(2) * 0.7, txt, 'Rotation', 90, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'Color', color, 'Interpreter', 'none');
end

function v = seriesRow(col, i)
%SERIESROW  Row i of an identifyMaterial series column, cell or numeric.
%   cell2table keeps orders/twoTheta as a cell column only when row vector
%   lengths differ; same-length rows (e.g. a single series) collapse to a
%   plain numeric matrix, where brace indexing would error.
    if iscell(col), v = col{i}; else, v = col(i, :); end
end

function s = shortName(name)
%SHORTNAME  Compact display alias for plot labels.
    map = struct('SrTiO3', "STO", 'SrRuO3', "SRO", 'PbTiO3', "PTO");
    n = char(name);
    if isfield(map, n), s = map.(n); else, s = string(name); end
end

function lines = identificationReport(R)
%IDENTIFICATIONREPORT  Results-panel text for an identifyMaterial run.
    lines = {};
    if R.substrate.found
        lines{end+1} = sprintf('Substrate %s: confirmed (%d peaks), c = %.4f A', ...
            R.substrate.name, numel(R.substrate.twoTheta), R.substrate.cMeas);
    else
        lines{end+1} = sprintf('Substrate %s: NOT FOUND in scan', R.substrate.name);
    end
    if ~isempty(R.ghosts) && height(R.ghosts) > 0
        lines{end+1} = sprintf('Filtered %d Kbeta/W-La ghost peak(s)', height(R.ghosts));
    end
    for i = 1:height(R.series)
        lines{end+1} = ''; %#ok<AGROW>
        lines{end+1} = sprintf('Series %d: c = %.4f A (orders: %s)', ...
            i, R.series.cMeas(i), strjoin(string(seriesRow(R.series.orders, i)), ',')); %#ok<AGROW>
        cand = R.series.candidates{i};
        if isempty(cand)
            lines{end+1} = '  -> no database match (unidentified)'; %#ok<AGROW>
        end
        for k = 1:numel(cand)
            extra = '';
            if ~isnan(cand(k).x)
                extra = sprintf('   x(Zr) ~ %.2f', cand(k).x);
            end
            lines{end+1} = sprintf( ...
                '  -> %-7s score %.2f   strain vs bulk %+.2f%%   relax %.2f%s', ...
                cand(k).name, cand(k).score, 100 * cand(k).strainVsBulk, ...
                cand(k).relaxation, extra); %#ok<AGROW>
        end
        fl = R.series.flags{i};
        if ~isempty(fl)
            lines{end+1} = sprintf('  flags: %s', strjoin(fl, ', ')); %#ok<AGROW>
        end
    end
    if ~isempty(R.unassigned)
        lines{end+1} = '';
        lines{end+1} = sprintf('Unassigned peaks: %s', ...
            strjoin(compose('%.2f', R.unassigned), ', '));
    end
    for nt = R.notes(:).'
        lines{end+1} = ''; %#ok<AGROW>
        lines{end+1} = char("NOTE: " + nt); %#ok<AGROW>
    end
end

function runXRR(fig)
    st   = fig.UserData;
    scan = st.scan; ax = st.ax;

    semilogy(ax, scan.twoTheta, max(scan.counts, 1), '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5);
    stylePubAxes(ax, '2\theta (°)', 'Counts', 'XRR');
    xlim(ax, [0, min(5, scan.twoTheta(end))]);

    % Full XRR pipeline (critical-edge → polynomial detrend → FFT period
    % → fringe peak detection → sinθ linear fit) lives in
    % xrdc.xrr.analyzeFringes. See that file for the algorithm; this
    % wrapper just exposes the GUI parameter panel.
    xrrMin  = getNum(st.params, 'xrrMin',   0);     % 0 = auto critical edge
    xrrMax  = getNum(st.params, 'xrrMax',   5.0);
    promPct = getNum(st.params, 'xrrProm',  1.5);   % % ripple on log decade

    args = {'UpperBound', xrrMax, 'MinProminence', promPct / 100};
    if xrrMin > 0
        args = [args, {'LowerBound', xrrMin}];
    end

    try
        res = xrdc.xrr.analyzeFringes(scan, args{:});
    catch ME
        writeResults(fig, {sprintf('XRR analysis failed: %s', ME.message)});
        return
    end

    if xrrMin > 0
        edgeMsg = sprintf('Lower bound (manual) = %.3f°', res.lowerBound);
    else
        edgeMsg = sprintf('Critical edge (auto) ≈ %.3f° → start at %.3f°', ...
            res.twoThetaC, res.lowerBound);
    end
    lines = {edgeMsg, ...
        sprintf('Fringe search range: [%.3f°, %.3f°]', ...
            res.lowerBound, res.upperBound), ...
        sprintf('Fringes detected: %d', res.nFringesDetected)};

    if ~isnan(res.thicknessFftNm)
        lines{end+1} = sprintf('FFT period: %.3f° (SNR %.1f)', ...
            res.fringePeriodDeg, res.fringeSnr);
    end

    if res.nFringesDetected >= 3
        hold(ax, 'on');
        plot(ax, [res.peaks.twoTheta], [res.peaks.counts], 'v', ...
            'MarkerFaceColor', [0.85 0.2 0.2], 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
        hold(ax, 'off');
        title(ax, sprintf('XRR — d = %.1f ± %.1f nm', ...
            res.thicknessQuadNm, res.thicknessQuadSeNm));

        lines = [lines, { ...
            '', ...
            sprintf('Thickness (Kiessig quad) = %.2f ± %.2f nm', ...
                res.thicknessQuadNm, res.thicknessQuadSeNm), ...
            sprintf('Thickness (linear sinθ)  = %.2f ± %.2f nm', ...
                res.thicknessFitNm, res.thicknessFitSeNm), ...
            sprintf('Thickness (FFT)          = %.2f nm', res.thicknessFftNm), ...
            sprintf('Recovered θ_c            = %.3f° (TER plateau %.3f°)', ...
                res.twoThetaCRecovered, res.twoThetaC), ...
            sprintf('λ                        = %.4f Å', scan.lambda)}];
    elseif ~isnan(res.thicknessFftNm)
        title(ax, sprintf('XRR — d ≈ %.1f nm (FFT only)', res.thicknessFftNm));
        lines = [lines, { ...
            '', ...
            sprintf('Thickness (FFT) = %.2f nm', res.thicknessFftNm), ...
            '— need ≥3 detected fringes for the quadratic Kiessig fit', ...
            sprintf('λ               = %.4f Å', scan.lambda)}];
    else
        lines{end+1} = '— could not determine thickness';
    end
    writeResults(fig, lines);
end

function runPhiScan(fig)
    st   = fig.UserData;
    scan = st.scan; ax = st.ax;

    % phi-aware pole detection: a NOISE-floor threshold (robust on weak
    % films, where "10% of max" sat below the noise and counted spikes as
    % poles) plus circular handling of the 360° wrap. See
    % xrdc.peaks.findPhiPeaks.
    nSig = getNum(st.params, 'noiseSigmas', 6);
    [pk, info] = xrdc.peaks.findPhiPeaks(scan, 'NoiseSigmas', nSig);

    plot(ax, scan.twoTheta, scan.counts, '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5);
    hold(ax, 'on');
    if ~isempty(pk)
        plot(ax, [pk.twoTheta], [pk.counts], 'v', ...
            'MarkerFaceColor', [0.85 0.2 0.2], 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    % Wrap repeats: same crystal direction as a counted pole, 360° away (the
    % cut-off edge of a pole whose apex is just past the scan). Mark with a
    % faint open grey marker + "↻ wrap" so a viewer sees the full turn was
    % covered, clearly distinct from the solid red pole markers — and NOT
    % counted toward the pole/symmetry tally.
    wr = info.wrapRepeats;
    if ~isempty(wr)
        plot(ax, [wr.twoTheta], [wr.counts], 'o', ...
            'MarkerEdgeColor', [0.55 0.55 0.55], 'MarkerFaceColor', 'none', ...
            'LineWidth', 1.2, 'MarkerSize', 11);
        for j = 1:numel(wr)
            text(ax, wr(j).twoTheta, double(wr(j).counts), '↻ wrap  ', ...
                'Color', [0.55 0.55 0.55], 'FontSize', 9, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
        end
    end
    hold(ax, 'off');
    stylePubAxes(ax, '\phi (°)', 'Counts', sprintf('φ scan — %d pole(s)', numel(pk)));
    set(ax, 'YScale', 'linear');

    lines = {sprintf('Detected %d pole(s)  (height ≥ %.0f cts; bg %.1f, noise σ %.1f):', ...
        numel(pk), info.heightThreshold, info.background, info.noiseSigma), ''};
    for i = 1:numel(pk)
        lines{end+1} = sprintf('  %2d.  φ = %7.2f°    I = %.0f', ...
            i, pk(i).twoTheta, pk(i).counts); %#ok<AGROW>
    end
    if ~isnan(info.fold)
        lines{end+1} = '';
        lines{end+1} = sprintf('Unique poles: %d    spacings: %s°', ...
            info.nUnique, join(string(round(info.spacings, 1)), ', '));
        lines{end+1} = sprintf('→ %d-fold symmetry', info.fold);
    end
    if ~isempty(wr)
        lines{end+1} = '';
        lines{end+1} = sprintf(['%d wrap repeat(s) at %s° (= a counted pole, ' ...
            '360° away; not counted)'], numel(wr), ...
            join(string(round([wr.twoTheta], 1)), ', '));
    end
    writeResults(fig, lines);
end

function runRSM(fig)
    st    = fig.UserData;
    ax    = st.ax;
    scans = st.rsmScans;

    if isempty(scans)
        writeResults(fig, {'No RSM slices loaded.'}); return
    end

    imin   = getNum(st.params, 'imin',      10);
    imax   = getNum(st.params, 'imax',      1e5);
    nCont  = getNum(st.params, 'nContours', 40);
    cmap   = getStr(st.params, 'colormap',  'turbo');

    xrdc.plot.plotRsm(scans, ...
        'TargetAxes', ax, ...
        'Mode',       "contourf", ...
        'NContours',  nCont, ...
        'Imin',       imin, ...
        'Imax',       imax, ...
        'Colormap',   string(cmap));

    writeResults(fig, { ...
        sprintf('Slices       : %d', numel(scans)), ...
        sprintf('ω range      : [%.3f°, %.3f°]', min([scans.secondAxis]), max([scans.secondAxis])), ...
        sprintf('2θ range     : [%.3f°, %.3f°]', scans(1).twoTheta(1), scans(1).twoTheta(end)), ...
        sprintf('Points/slice : %d', numel(scans(1).twoTheta)), ...
        '', ...
        sprintf('Colorbar     : [%g, %g]', imin, imax), ...
        sprintf('Colormap     : %s', cmap)});
end

function plotBasic(fig)
    st   = fig.UserData;
    scan = st.scan; ax = st.ax;
    semilogy(ax, scan.twoTheta, max(scan.counts, 1), '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5);
    stylePubAxes(ax, '2\theta (°)', 'Counts', char(scan.identifier));
    writeResults(fig, { ...
        sprintf('Identifier : %s', scan.identifier), ...
        sprintf('Scan type  : %s', scan.scanType), ...
        sprintf('Points     : %d', numel(scan.twoTheta)), ...
        sprintf('Wavelength : %.4f Å', scan.lambda)});
end

% =====================================================================
% Helpers
% =====================================================================
function stylePubAxes(ax, xlab, ylab, ttl)
    xlabel(ax, xlab); ylabel(ax, ylab);
    if ~isempty(ttl), title(ax, ttl); end
    grid(ax, 'on');
    set(ax, 'FontName', 'Arial', 'FontSize', 13, 'LineWidth', 1.25, 'Box', 'on');
    ax.Title.FontSize = 15;
    ax.XLabel.FontSize = 14; ax.YLabel.FontSize = 14;
end

function writeResults(fig, lines)
    st = fig.UserData;
    if ~isempty(st.resultsArea) && isvalid(st.resultsArea)
        st.resultsArea.Value = lines(:);
    end
end

% =====================================================================
% Plot-style customization (publication export)
% =====================================================================
function s = defaultStyle()
%DEFAULTSTYLE  Override fields, all empty/"auto" → keep analysis defaults.
    s = struct('title', '', 'xlabel', '', 'ylabel', '', ...
        'xmin', '', 'xmax', '', 'ymin', '', 'ymax', '', ...
        'yscale', 'auto', 'fontSize', '', 'lineWidth', '', ...
        'lineColor', 'auto', 'markers', 'auto', 'grid', 'auto', ...
        'exportW', '', 'exportH', '');
end

function v = styleNum(style, name)
    if isfield(style, name)
        v = str2double(strtrim(string(style.(name))));
    else
        v = NaN;
    end
end

function onCustomizePlot(fig)
%ONCUSTOMIZEPLOT  Modal-ish dialog of optional overrides for the live plot.
%   Defaults are preserved: blank/"auto" fields change nothing. "Apply"
%   re-renders the current analysis with the overrides; "Reset" restores
%   defaults. Export honours the width/height fields.
    st = fig.UserData;
    if isempty(st.scan) && isempty(st.rsmScans), return, end
    s = st.style;
    T = appTheme();

    d = uifigure('Name', 'Customize Plot', 'Position', [220 160 380 590], ...
        'Color', T.bg);
    try, d.Theme = 'light'; catch, end
    g = uigridlayout(d, [17 2]);
    g.RowHeight   = [repmat({28}, 1, 16), {36}];
    g.ColumnWidth = {150, '1x'};
    g.Padding     = [12 12 12 12];
    g.RowSpacing  = 5;
    g.BackgroundColor = T.bg;

    f = struct(); r = 1;
    [f.title, r]     = dlgEdit(g, r, 'Title (blank=auto)', s.title);
    [f.xlabel, r]    = dlgEdit(g, r, 'X label',            s.xlabel);
    [f.ylabel, r]    = dlgEdit(g, r, 'Y label',            s.ylabel);
    [f.xmin, r]      = dlgEdit(g, r, 'X min',              s.xmin);
    [f.xmax, r]      = dlgEdit(g, r, 'X max',              s.xmax);
    [f.ymin, r]      = dlgEdit(g, r, 'Y min',              s.ymin);
    [f.ymax, r]      = dlgEdit(g, r, 'Y max',              s.ymax);
    [f.yscale, r]    = dlgDrop(g, r, 'Y scale', {'auto','linear','log'}, s.yscale);
    [f.fontSize, r]  = dlgEdit(g, r, 'Font size',          s.fontSize);
    [f.lineWidth, r] = dlgEdit(g, r, 'Line width',         s.lineWidth);
    [f.lineColor, r] = dlgDrop(g, r, 'Line colour', ...
        {'auto','blue','black','red','green','orange','purple','gray'}, s.lineColor);
    [f.markers, r]   = dlgDrop(g, r, 'Peak markers', {'auto','on','off'}, s.markers);
    [f.grid, r]      = dlgDrop(g, r, 'Grid',         {'auto','on','off'}, s.grid);
    [f.exportW, r]   = dlgEdit(g, r, 'Export width (in)',  s.exportW);
    [f.exportH, r]   = dlgEdit(g, r, 'Export height (in)', s.exportH);

    note = uilabel(g, 'Text', 'Blank / "auto" keeps the default.', ...
        'FontAngle', 'italic', 'FontName', T.font, 'FontColor', T.textDim);
    note.Layout.Row = 16; note.Layout.Column = [1 2];

    btns = uigridlayout(g, [1 3], 'Padding', [0 0 0 0]);
    btns.Layout.Row = 17; btns.Layout.Column = [1 2];
    btns.BackgroundColor = T.bg;
    uibutton(btns, 'Text', 'Apply', 'FontWeight', 'bold', 'FontName', T.font, ...
        'BackgroundColor', T.gold, 'FontColor', T.ink, ...
        'ButtonPushedFcn', @(~,~) applyDlg());
    uibutton(btns, 'Text', 'Reset', 'FontName', T.font, ...
        'BackgroundColor', T.btn, 'FontColor', T.text, ...
        'ButtonPushedFcn', @(~,~) resetDlg());
    uibutton(btns, 'Text', 'Close', 'FontName', T.font, ...
        'BackgroundColor', T.btn, 'FontColor', T.text, ...
        'ButtonPushedFcn', @(~,~) close(d));

    function applyDlg()
        ns = defaultStyle();
        for fn = fieldnames(f)'
            ns.(fn{1}) = f.(fn{1}).Value;
        end
        s2 = fig.UserData; s2.style = ns; fig.UserData = s2;
        runAnalysis(fig);
    end
    function resetDlg()
        ds = defaultStyle();
        for fn = fieldnames(f)'
            f.(fn{1}).Value = ds.(fn{1});
        end
        s2 = fig.UserData; s2.style = ds; fig.UserData = s2;
        runAnalysis(fig);
    end
end

function [h, row] = dlgEdit(g, row, label, val)
    T = appTheme();
    lbl = uilabel(g, 'Text', label, 'FontName', T.font, 'FontColor', T.text);
    lbl.Layout.Row = row; lbl.Layout.Column = 1;
    h = uieditfield(g, 'Value', char(string(val)));
    h.Layout.Row = row; h.Layout.Column = 2;
    row = row + 1;
end

function [h, row] = dlgDrop(g, row, label, items, val)
    T = appTheme();
    lbl = uilabel(g, 'Text', label, 'FontName', T.font, 'FontColor', T.text);
    lbl.Layout.Row = row; lbl.Layout.Column = 1;
    h = uidropdown(g, 'Items', items, 'Value', char(string(val)));
    h.Layout.Row = row; h.Layout.Column = 2;
    row = row + 1;
end

function v = getNum(params, name, default)
    if isfield(params, name)
        v = str2double(params.(name));
        if isnan(v), v = default; end
    else
        v = default;
    end
end

function v = getStr(params, name, default)
    if isfield(params, name)
        v = char(params.(name));
    else
        v = default;
    end
end

function prefix = prefixForType(t)
    switch t
        case "omega",          prefix = "RC";
        case "twothetaomega",  prefix = "th2th";
        case "xrr",            prefix = "xrr";
        case "phi",            prefix = "phi";
        case "rsm",            prefix = "rsm";
        otherwise,             prefix = "scan";
    end
end

function placeholder(ax)
    cla(ax);
    text(ax, 0.5, 0.5, 'No scan loaded', ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 16, 'Color', [0.5 0.5 0.5]);
    set(ax, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none');
end
