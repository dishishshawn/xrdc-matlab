function publicationStyle(ax, options)
%PUBLICATIONSTYLE  Apply XRDC publication-ready styling to an axes.
%   xrdc.plot.publicationStyle()                        % current axes
%   xrdc.plot.publicationStyle(ax)
%   xrdc.plot.publicationStyle(ax, Name, Value, ...)
%
%   Implements the invariants from ALGORITHM_SPEC §10:
%     * Arial / Helvetica font
%     * 18 pt tick labels, 20 pt axis labels, 22 pt title
%     * line width 1.5 for lines already on the axes
%     * white background, minor ticks on, grid off
%     * XLabel `2θ (°)`, YLabel `Counts` (log-scale) by default
%
%   Style decisions are centralised here so every plotting function in
%   +xrdc/+plot/ stays visually consistent. Dr. Paik's paper-specific
%   overrides (font family, palette, tick sizes) can be applied by
%   passing Name/Value pairs.
%
%   Name/Value options
%     'FontName'       default "Arial"
%     'TickFontSize'   default 18
%     'LabelFontSize'  default 20
%     'TitleFontSize'  default 22
%     'LineWidth'      default 1.5 — applies to pre-existing Line children
%     'LogY'           default true
%     'Grid'           default false
%     'MinorTicks'     default true
%
%   See ALGORITHM_SPEC §10 for rationale and defaults.

    arguments
        ax                     (1,1) = gca
        options.FontName       (1,1) string  = "Arial"
        options.TickFontSize   (1,1) double  = 18
        options.LabelFontSize  (1,1) double  = 20
        options.TitleFontSize  (1,1) double  = 22
        options.LineWidth      (1,1) double  = 1.5
        options.LogY           (1,1) logical = true
        options.Grid           (1,1) logical = false
        options.MinorTicks     (1,1) logical = true
    end

    if ~isa(ax, 'matlab.graphics.axis.Axes')
        error('xrdc:plot:badAxes', ...
            'First argument must be a MATLAB Axes (use gca or a handle).');
    end

    set(ax, ...
        'FontName',            char(options.FontName), ...
        'FontSize',            options.TickFontSize, ...
        'Box',                 'on', ...
        'TickDir',             'in', ...
        'LineWidth',           1.25, ...      % axis line, not data line
        'Color',               'w');
    ax.XLabel.FontSize = options.LabelFontSize;
    ax.YLabel.FontSize = options.LabelFontSize;
    ax.Title.FontSize  = options.TitleFontSize;
    ax.XLabel.FontName = char(options.FontName);
    ax.YLabel.FontName = char(options.FontName);
    ax.Title.FontName  = char(options.FontName);

    if isempty(get(ax.XLabel, 'String'))
        xlabel(ax, '2\theta (\circ)');
    end
    if isempty(get(ax.YLabel, 'String'))
        ylabel(ax, 'Counts');
    end

    if options.LogY
        set(ax, 'YScale', 'log');
    else
        set(ax, 'YScale', 'linear');
    end

    if options.Grid
        grid(ax, 'on');
    else
        grid(ax, 'off');
    end

    if options.MinorTicks
        ax.XAxis.MinorTick = 'on';
        ax.YAxis.MinorTick = 'on';
    end

    % Style any existing Line children (data traces).
    lines = findobj(ax, 'Type', 'line');
    for i = 1:numel(lines)
        if lines(i).LineWidth < options.LineWidth
            lines(i).LineWidth = options.LineWidth;
        end
    end

    % White figure background for export.
    fig = ancestor(ax, 'figure');
    if ~isempty(fig)
        set(fig, 'Color', 'w');
    end
end
