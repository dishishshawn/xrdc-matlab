function h = plotScan(scan, options)
%PLOTSCAN  Publication-quality plot of a single XRD scan.
%   h = xrdc.plot.plotScan(scan)
%   h = xrdc.plot.plotScan(scan, Name, Value, ...)
%
%   Produces a figure that follows the invariants in ALGORITHM_SPEC §10
%   (log-Y, Arial 18 pt ticks, white background, etc.). The caller can
%   override any of those via the Name/Value pairs forwarded to
%   `xrdc.plot.publicationStyle`.
%
%   Behaviour:
%     - Clamps `counts <= 0` to 1 before plotting when LogY = true —
%       matches the Delphi PaintGraph clamp at xrdc1.pas:1901-1904 so
%       zeros do not send `log(0)` to -Inf.
%     - Optionally overlays detected peak markers if `scan.peaks` is a
%       non-empty struct array with a `.twoTheta` field.
%
%   Input
%     scan : struct from xrdc.io reader (needs .twoTheta, .counts).
%
%   Name/Value options (forwarded to publicationStyle unless noted)
%     'Color'          default [0 0.447 0.741] (MATLAB default blue)
%     'LineWidth'      default 1.5
%     'LogY'           default true
%     'Title'          default "" (falls back to scan.identifier)
%     'ShowPeaks'      default true if scan.peaks exists; marker set
%                      at (peak.twoTheta, peak.counts) with a down-
%                      triangle.
%     'TargetAxes'     default gca (pass a handle to draw into an
%                      existing figure, e.g. for a subplot).
%
%   Returns: struct with the line handle (.line), peak marker handle
%   (.peakMarkers), and axes handle (.ax) so callers can tweak further.
%
%   Note — final style details (colour palette, font size, export
%   resolution) may change once Dr. Paik provides a reference paper;
%   this function is structured so that changes flow through
%   `publicationStyle` without touching callers.

    arguments
        scan                     (1,1) struct
        options.Color            (1,3) double = [0 0.447 0.741]
        options.LineWidth        (1,1) double = 1.5
        options.LogY             (1,1) logical = true
        options.Title            (1,1) string = ""
        options.ShowPeaks        (1,1) logical = true
        options.TargetAxes                    = []
        % forwarded to publicationStyle
        options.FontName         (1,1) string  = "Arial"
        options.TickFontSize     (1,1) double  = 18
        options.LabelFontSize    (1,1) double  = 20
        options.TitleFontSize    (1,1) double  = 22
    end

    if ~isfield(scan, 'twoTheta') || ~isfield(scan, 'counts')
        error('xrdc:plot:badScan', ...
            'scan must have .twoTheta and .counts fields.');
    end

    if isempty(options.TargetAxes)
        fig = figure();
        ax  = axes(fig);
    else
        ax = options.TargetAxes;
        fig = ancestor(ax, 'figure');
    end

    x = double(scan.twoTheta(:));
    y = double(scan.counts(:));
    if options.LogY
        % Delphi clamp: counts <= 0 → 1 (xrdc1.pas:1901-1904).
        y(y <= 0) = 1;
    end

    lineHandle = plot(ax, x, y, ...
        'Color',     options.Color, ...
        'LineWidth', options.LineWidth);
    hold(ax, 'on');

    % Overlay peak markers if present
    peakMarkers = [];
    if options.ShowPeaks && isfield(scan, 'peaks') && ~isempty(scan.peaks)
        pk = scan.peaks;
        if isfield(pk, 'twoTheta') && isfield(pk, 'counts')
            xp = [pk.twoTheta];
            yp = [pk.counts];
            if options.LogY
                yp(yp <= 0) = 1;
            end
            peakMarkers = plot(ax, xp, yp, 'v', ...
                'MarkerEdgeColor', [0.6 0 0], ...
                'MarkerFaceColor', [0.85 0.35 0.35], ...
                'MarkerSize',      8, ...
                'LineStyle',       'none');
        end
    end

    % Title: explicit > scan.identifier > nothing
    if strlength(options.Title) > 0
        title(ax, options.Title);
    elseif isfield(scan, 'identifier') && strlength(scan.identifier) > 0
        title(ax, char(scan.identifier));
    end
    hold(ax, 'off');

    xrdc.plot.publicationStyle(ax, ...
        'FontName',       options.FontName, ...
        'TickFontSize',   options.TickFontSize, ...
        'LabelFontSize',  options.LabelFontSize, ...
        'TitleFontSize',  options.TitleFontSize, ...
        'LineWidth',      options.LineWidth, ...
        'LogY',           options.LogY);

    h = struct( ...
        'line',        lineHandle, ...
        'peakMarkers', peakMarkers, ...
        'ax',          ax, ...
        'figure',      fig);
end
