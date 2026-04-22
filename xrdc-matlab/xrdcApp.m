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

    fig = uifigure('Name', 'XRDC Scan Analyzer', 'Position', [100 100 1200 750]);
    grid = uigridlayout(fig, [3 2]);
    grid.RowHeight    = {40, 32, '1x'};
    grid.ColumnWidth  = {300, '1x'};
    grid.RowSpacing   = 6;
    grid.ColumnSpacing = 6;
    grid.Padding      = [8 8 8 8];

    % Top bar: load + export
    topBar = uigridlayout(grid, [1 3]);
    topBar.Layout.Row = 1; topBar.Layout.Column = [1 2];
    topBar.ColumnWidth = {130, 130, '1x'};
    topBar.Padding = [0 0 0 0];
    uibutton(topBar, 'Text', 'Load Scan...', ...
        'FontSize', 13, 'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(~,~) onLoadScan(fig));
    exportBtn = uibutton(topBar, 'Text', 'Export 600 dpi...', ...
        'FontSize', 13, ...
        'Enable', 'off', ...
        'ButtonPushedFcn', @(~,~) onExport(fig));

    % Info strip
    infoLbl = uilabel(grid, 'Text', '  No scan loaded. Click "Load Scan..." to begin.', ...
        'FontSize', 12, 'FontColor', [0.4 0.4 0.4], ...
        'HorizontalAlignment', 'left');
    infoLbl.Layout.Row = 2; infoLbl.Layout.Column = [1 2];

    % Left: analysis panel
    leftPanel = uipanel(grid, 'Title', 'Analysis', 'FontSize', 13, 'FontWeight', 'bold');
    leftPanel.Layout.Row = 3; leftPanel.Layout.Column = 1;

    % Right: plot
    plotPanel = uipanel(grid, 'Title', 'Preview', 'FontSize', 13, 'FontWeight', 'bold');
    plotPanel.Layout.Row = 3; plotPanel.Layout.Column = 2;
    plotGrid = uigridlayout(plotPanel, [1 1]);
    plotGrid.Padding = [4 4 4 4];
    ax = uiaxes(plotGrid);

    % Store state on the figure so callbacks can share it
    st = struct();
    st.scan        = [];
    st.rsmScans    = [];
    st.filePath    = '';
    st.detectedType = "";
    st.params      = struct();
    st.ax          = ax;
    st.infoLbl     = infoLbl;
    st.exportBtn   = exportBtn;
    st.leftPanel   = leftPanel;
    st.resultsArea = [];
    fig.UserData = st;

    placeholder(ax);
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
    st.rsmScans = [];

    try
        isRsm = endsWith(lower(file), '.xrdml') && contains(lower(file), 'rsm');
        if isRsm
            st.rsmScans    = xrdc.rsm.loadAreaScan({fullPath});
            st.scan        = st.rsmScans(1);   % representative
            st.detectedType = "rsm";
        else
            st.scan        = xrdc.io.readScan(fullPath);
            st.detectedType = lower(string(st.scan.scanType));
            % XRR is stored as twoThetaOmega; disambiguate by filename
            if st.detectedType == "twothetaomega" && contains(lower(file), "xrr")
                st.detectedType = "xrr";
            end
        end
    catch ME
        uialert(fig, sprintf('Failed to load the file:\n\n%s', ME.message), ...
            'Load error', 'Icon', 'error');
        return
    end

    st.infoLbl.Text = sprintf('  %s     |     %s     |     %d points', ...
        file, upper(st.detectedType), numel(st.scan.twoTheta));
    st.infoLbl.FontColor = [0 0 0];
    st.exportBtn.Enable = 'on';
    fig.UserData = st;

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

    try
        exportgraphics(st.ax, fullfile(path, file), 'Resolution', 600);
        uialert(fig, sprintf('Saved:\n%s', fullfile(path, file)), ...
            'Export', 'Icon', 'success');
    catch ME
        uialert(fig, sprintf('Export failed:\n%s', ME.message), 'Error', 'Icon', 'error');
    end
end

function onParamChange(fig, name, value)
    st = fig.UserData;
    st.params.(name) = value;
    fig.UserData = st;
    runAnalysis(fig);
end

% =====================================================================
% Panel construction per scan type
% =====================================================================
function buildAnalysisPanel(fig)
    st = fig.UserData;
    delete(st.leftPanel.Children);

    g = uigridlayout(st.leftPanel, [14 2]);
    g.RowHeight   = [repmat({26}, 1, 6), {'1x'}, repmat({26}, 1, 7)];
    g.ColumnWidth = {110, '1x'};
    g.RowSpacing  = 4;
    g.Padding     = [6 6 6 6];

    row = 1;
    switch st.detectedType
        case "omega"
            row = addEdit (g, row, 'Fit window (°)', '0.5',      @(v) onParamChange(fig, 'fitWindow', v));
            row = addDrop (g, row, 'Shape',          {'lorentz','gauss','pseudoVoigt'}, 'lorentz', ...
                                                                   @(v) onParamChange(fig, 'shape',     v));
        case "twothetaomega"
            row = addEdit (g, row, 'Min prom (%)',   '5',        @(v) onParamChange(fig, 'promPct',   v));
        case "xrr"
            row = addEdit (g, row, 'Fringe 2θ min',  '0.5',      @(v) onParamChange(fig, 'xrrMin',    v));
            row = addEdit (g, row, 'Fringe 2θ max',  '3.0',      @(v) onParamChange(fig, 'xrrMax',    v));
            row = addEdit (g, row, 'Min prom (%)',   '2',        @(v) onParamChange(fig, 'xrrProm',   v));
        case "phi"
            row = addEdit (g, row, 'Min prom (%)',   '10',       @(v) onParamChange(fig, 'promPct',   v));
        case "rsm"
            row = addDrop (g, row, 'Colormap',       {'turbo','parula','jet'}, 'turbo', ...
                                                                   @(v) onParamChange(fig, 'colormap', v));
            row = addEdit (g, row, 'Imin (counts)',  '1',        @(v) onParamChange(fig, 'imin',      v));
            row = addEdit (g, row, 'Imax (counts)',  '1e5',      @(v) onParamChange(fig, 'imax',      v));
            row = addEdit (g, row, 'Contours',       '40',       @(v) onParamChange(fig, 'nContours', v));
    end

    % Results area fills the rest
    hdr = uilabel(g, 'Text', 'Results', 'FontWeight', 'bold');
    hdr.Layout.Row = 8; hdr.Layout.Column = [1 2];
    ta = uitextarea(g, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 11);
    ta.Layout.Row = [9 14]; ta.Layout.Column = [1 2];
    st.resultsArea = ta;
    fig.UserData = st;
end

function row = addEdit(g, row, label, default, cb)
    lbl = uilabel(g, 'Text', label);
    lbl.Layout.Row = row; lbl.Layout.Column = 1;
    ef = uieditfield(g, 'Value', default, ...
        'ValueChangedFcn', @(src, ~) cb(src.Value));
    ef.Layout.Row = row; ef.Layout.Column = 2;
    row = row + 1;
end

function row = addDrop(g, row, label, items, default, cb)
    lbl = uilabel(g, 'Text', label);
    lbl.Layout.Row = row; lbl.Layout.Column = 1;
    dd = uidropdown(g, 'Items', items, 'Value', default, ...
        'ValueChangedFcn', @(src, ~) cb(src.Value));
    dd.Layout.Row = row; dd.Layout.Column = 2;
    row = row + 1;
end

% =====================================================================
% Analysis dispatchers
% =====================================================================
function runAnalysis(fig)
    st = fig.UserData;
    cla(st.ax); reset(st.ax);
    try
        switch st.detectedType
            case "omega",         runRockingCurve(fig);
            case "twothetaomega", runThetaTwoTheta(fig);
            case "xrr",           runXRR(fig);
            case "phi",           runPhiScan(fig);
            case "rsm",           runRSM(fig);
            otherwise,            plotBasic(fig);
        end
    catch ME
        uialert(fig, sprintf('Analysis error:\n\n%s', ME.message), ...
            'Analysis error', 'Icon', 'error');
    end
end

function runRockingCurve(fig)
    st   = fig.UserData;
    scan = st.scan;    ax = st.ax;

    pk = xrdc.peaks.findPeaks(scan, 'MinProminence', max(scan.counts) * 0.05);
    semilogy(ax, scan.twoTheta, max(scan.counts, 1), '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5);
    stylePubAxes(ax, '\omega (°)', 'Counts', '');

    if isempty(pk)
        writeResults(fig, {'No peak detected — try lowering threshold.'});
        return
    end

    [~, idx] = max([pk.counts]); pkMain = pk(idx);

    w     = getNum(st.params, 'fitWindow', 0.5);
    shape = getStr(st.params, 'shape',    'lorentz');
    window = [pkMain.twoTheta - w, pkMain.twoTheta + w];

    fit = xrdc.peaks.fitPeak(scan, window, 'Shape', string(shape));

    hold(ax, 'on');
    plot(ax, fit.xFit, max(fit.yFit, 1), '--', ...
        'Color', [0.85 0.2 0.2], 'LineWidth', 1.8);
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
        '', ...
        sprintf('Fit window: ±%.3f° around peak', w)});
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

function runXRR(fig)
    st   = fig.UserData;
    scan = st.scan; ax = st.ax;

    semilogy(ax, scan.twoTheta, max(scan.counts, 1), '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5);
    stylePubAxes(ax, '2\theta (°)', 'Counts', 'XRR');
    xlim(ax, [0, min(5, scan.twoTheta(end))]);

    xrrMin = getNum(st.params, 'xrrMin',   0.5);
    xrrMax = getNum(st.params, 'xrrMax',   3.0);
    prom   = getNum(st.params, 'xrrProm',  2);

    mask = scan.twoTheta > xrrMin & scan.twoTheta < xrrMax;
    sub  = scan; sub.twoTheta = scan.twoTheta(mask); sub.counts = scan.counts(mask);

    lines = {sprintf('Fringe search range: [%.2f°, %.2f°]', xrrMin, xrrMax)};
    if nnz(mask) < 5
        lines{end+1} = '— not enough points in range';
        writeResults(fig, lines); return
    end

    pk = xrdc.peaks.findPeaks(sub, ...
        'MinProminence', max(sub.counts) * prom / 100, ...
        'MinSeparation', 0.05);
    lines{end+1} = sprintf('Fringes detected: %d', numel(pk));

    if numel(pk) >= 2
        ttPk  = [pk.twoTheta];
        thick = xrdc.lattice.thicknessFromFringes(ttPk(:), scan.lambda);
        t_nm  = thick.thicknessFitNm;

        hold(ax, 'on');
        plot(ax, [pk.twoTheta], [pk.counts], 'v', ...
            'MarkerFaceColor', [0.85 0.2 0.2], 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
        hold(ax, 'off');
        title(ax, sprintf('XRR — d = %.1f ± %.1f nm', t_nm, thick.thicknessFitSeNm));

        lines = [lines, { ...
            '', ...
            sprintf('Thickness  = %.2f nm', t_nm), ...
            sprintf('Uncertainty = ± %.2f nm', thick.thicknessFitSeNm), ...
            sprintf('λ          = %.4f Å', scan.lambda)}];
    else
        lines{end+1} = '— need ≥2 fringes for thickness';
    end
    writeResults(fig, lines);
end

function runPhiScan(fig)
    st   = fig.UserData;
    scan = st.scan; ax = st.ax;

    promPct = getNum(st.params, 'promPct', 10);
    pk = xrdc.peaks.findPeaks(scan, ...
        'MinProminence', max(scan.counts) * promPct / 100, ...
        'MinSeparation', 30);

    plot(ax, scan.twoTheta, scan.counts, '-', 'Color', [0.1 0.4 0.8], 'LineWidth', 1.5);
    hold(ax, 'on');
    if ~isempty(pk)
        plot(ax, [pk.twoTheta], [pk.counts], 'v', ...
            'MarkerFaceColor', [0.85 0.2 0.2], 'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    hold(ax, 'off');
    stylePubAxes(ax, '\phi (°)', 'Counts', sprintf('φ scan — %d peaks', numel(pk)));
    set(ax, 'YScale', 'linear');

    lines = {sprintf('Detected %d peak(s):', numel(pk)), ''};
    for i = 1:numel(pk)
        lines{end+1} = sprintf('  %2d.  φ = %7.2f°    I = %.0f', ...
            i, pk(i).twoTheta, pk(i).counts); %#ok<AGROW>
    end
    if numel(pk) >= 4
        phis     = sort([pk.twoTheta]);
        spacings = diff(phis);
        lines{end+1} = '';
        lines{end+1} = sprintf('Spacings: %s °', join(string(round(spacings, 1)), ', '));
        if numel(pk) == 4 && all(abs(spacings - 90) < 5)
            lines{end+1} = '✓ 4-fold symmetry (within 5°)';
        end
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

    imin   = getNum(st.params, 'imin',      1);
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
