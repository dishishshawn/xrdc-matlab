function h = plotStack(scans, options)
%PLOTSTACK  Waterfall plot of multiple scans on a shared axis.
%   h = xrdc.plot.plotStack(scans)
%   h = xrdc.plot.plotStack(scans, Name, Value, ...)
%
%   Stacks multiple XRD scans vertically using the Delphi "multiplier"
%   convention (xrdc1.pas:1939-1953):
%     * For N ≤ 50 scans, scan j (0-indexed) is scaled by `3^(j+1)`.
%     * For N > 50 scans, scaled by `10^(20·(j+1)/N)` instead, so the
%       stack fits on a log-Y axis without running off the top.
%
%   The scaling is purely cosmetic — peak *positions* are unchanged,
%   but intensities on the y-axis are no longer meaningful after
%   stacking, so set `YLabel` to "Intensity (arb.)" downstream.
%
%   Input
%     scans : struct array (N×1) of scan structs. Each must have
%             .twoTheta and .counts; .identifier is used for the
%             legend if present.
%
%   Name/Value options
%     'Palette'    default "lines" — a MATLAB colormap name or an
%                  N×3 matrix of RGB rows. Replaces the Delphi random
%                  per-scan colour (see ALGORITHM_SPEC §10 rationale).
%     'LogY'       default true
%     'LineWidth'  default 1.2 (1.5 for the first trace)
%     'Title'      default ""
%     'Multiplier' override the auto-multiplier with an Nx1 vector or
%                  a scalar base factor. Pass [] to use the default.
%     'TargetAxes' default [] (new figure)
%     Additional Name/Value pairs forwarded to publicationStyle.
%
%   Returns a struct with .lines (Nx1 handle array), .ax, .figure.

    arguments
        scans                             struct
        options.Palette                   = "lines"
        options.LogY            (1,1) logical = true
        options.LineWidth       (1,1) double = 1.2
        options.Title           (1,1) string = ""
        options.Multiplier                  = []
        options.TargetAxes                  = []
        options.FontName        (1,1) string = "Arial"
        options.TickFontSize    (1,1) double = 18
        options.LabelFontSize   (1,1) double = 20
        options.TitleFontSize   (1,1) double = 22
    end

    N = numel(scans);
    if N == 0
        error('xrdc:plot:emptyStack', 'plotStack requires at least one scan.');
    end

    % Default multiplier follows the Delphi convention.
    if isempty(options.Multiplier)
        if N <= 50
            mult = 3 .^ (0:N-1);          % first scan unscaled
        else
            mult = 10 .^ (20 * (0:N-1) / N);
        end
    elseif isscalar(options.Multiplier)
        mult = options.Multiplier .^ (0:N-1);
    else
        mult = options.Multiplier(:).';
        if numel(mult) ~= N
            error('xrdc:plot:badMultiplier', ...
                'Multiplier vector must have %d elements (got %d).', N, numel(mult));
        end
    end

    % Resolve palette to an Nx3 matrix.
    cmap = resolvePalette(options.Palette, N);

    if isempty(options.TargetAxes)
        fig = figure();
        ax  = axes(fig);
    else
        ax  = options.TargetAxes;
        fig = ancestor(ax, 'figure');
    end

    lineHandles = gobjects(N, 1);
    legendTexts = cell(N, 1);
    hold(ax, 'on');
    for j = 1:N
        s = scans(j);
        x = double(s.twoTheta(:));
        y = double(s.counts(:)) * mult(j);
        if options.LogY
            y(y <= 0) = 1;
        end
        lw = options.LineWidth;
        if j == 1, lw = max(lw, 1.5); end
        lineHandles(j) = plot(ax, x, y, ...
            'Color',     cmap(j, :), ...
            'LineWidth', lw);
        if isfield(s, 'identifier') && strlength(s.identifier) > 0
            legendTexts{j} = char(s.identifier);
        else
            legendTexts{j} = sprintf('Scan %d', j);
        end
    end
    hold(ax, 'off');
    if strlength(options.Title) > 0
        title(ax, options.Title);
    end
    legend(ax, legendTexts, 'Location', 'best');

    xrdc.plot.publicationStyle(ax, ...
        'FontName',       options.FontName, ...
        'TickFontSize',   options.TickFontSize, ...
        'LabelFontSize',  options.LabelFontSize, ...
        'TitleFontSize',  options.TitleFontSize, ...
        'LineWidth',      options.LineWidth, ...
        'LogY',           options.LogY);

    % Intensity axis is no longer meaningful in absolute terms.
    ylabel(ax, 'Intensity (arb.)');

    h = struct('lines', lineHandles, 'ax', ax, 'figure', fig);
end

% -------------------------------------------------------------------------

function cmap = resolvePalette(pal, N)
    if isnumeric(pal)
        if size(pal, 2) ~= 3
            error('xrdc:plot:badPalette', ...
                'Numeric palette must have 3 columns (RGB).');
        end
        if size(pal, 1) < N
            % Tile if too few colours
            reps = ceil(N / size(pal, 1));
            cmap = repmat(pal, reps, 1);
        else
            cmap = pal;
        end
    else
        name = char(pal);
        try
            fn = str2func(name);
            cmap = fn(N);
        catch
            error('xrdc:plot:badPalette', ...
                'Unknown palette "%s" (expected a colormap name or Nx3 matrix).', name);
        end
    end
    cmap = cmap(1:N, :);
end
