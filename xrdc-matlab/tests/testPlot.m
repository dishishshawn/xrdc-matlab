function tests = testPlot
%TESTPLOT  Unit tests for the +xrdc.+plot package.
%   These tests exercise the plotting API — they create figures,
%   assert on returned handle structs, and close the figure in
%   teardown. Visual verification of the publication style is the
%   user's responsibility (and will be refined once Dr. Paik picks
%   a reference paper).
    tests = functiontests(localfunctions);
end

% ---------- fixtures ----------

function setupOnce(testCase) %#ok<INUSD>
    % Use a CurrentFigureFixture-equivalent: don't steal focus
    set(groot, 'DefaultFigureVisible', 'off');
end

function teardownOnce(testCase) %#ok<INUSD>
    set(groot, 'DefaultFigureVisible', 'on');
    close all force
end

function s = makeSyntheticScan()
    x = (20:0.02:80).';
    centers = [28, 47, 69];
    fwhms   = [0.3, 0.25, 0.35];
    amps    = [8000, 6000, 10000];
    y = 150 * ones(size(x));
    for k = 1:numel(centers)
        sigma = fwhms(k) / (2 * sqrt(2 * log(2)));
        y = y + amps(k) * exp(-(x - centers(k)).^2 ./ (2 * sigma^2));
    end
    s = xrdc.io.emptyScan();
    s.twoTheta   = x;
    s.counts     = y;
    s.identifier = "synthetic";
end

% ---------- plotScan ----------

function testPlotScanReturnsHandles(tc)
    s = makeSyntheticScan();
    h = xrdc.plot.plotScan(s);
    tc.verifyTrue(isa(h.ax, 'matlab.graphics.axis.Axes'));
    tc.verifyTrue(isa(h.line, 'matlab.graphics.chart.primitive.Line'));
    tc.verifyTrue(isgraphics(h.figure, 'figure'));
    close(h.figure);
end

function testPlotScanLogYClampsZeros(tc)
    s = makeSyntheticScan();
    s.counts(1:5) = 0;                           % inject zeros
    h = xrdc.plot.plotScan(s, 'LogY', true);
    ydata = h.line.YData;
    tc.verifyTrue(all(ydata > 0), ...
        'LogY=true should clamp non-positive counts to 1.');
    close(h.figure);
end

function testPlotScanLinearYAllowsZeros(tc)
    s = makeSyntheticScan();
    s.counts(1:5) = 0;
    h = xrdc.plot.plotScan(s, 'LogY', false);
    ydata = h.line.YData;
    tc.verifyEqual(ydata(1:5), zeros(1,5), 'AbsTol', 0);
    close(h.figure);
end

function testPlotScanUsesIdentifierAsTitle(tc)
    s = makeSyntheticScan();
    s.identifier = "SrTiO3 (001)";
    h = xrdc.plot.plotScan(s);
    tc.verifyEqual(string(h.ax.Title.String), "SrTiO3 (001)");
    close(h.figure);
end

function testPlotScanOverlaysPeaks(tc)
    s = makeSyntheticScan();
    s.peaks = [struct('twoTheta', 28, 'counts', 8150); ...
               struct('twoTheta', 47, 'counts', 6150); ...
               struct('twoTheta', 69, 'counts', 10150)];
    h = xrdc.plot.plotScan(s, 'ShowPeaks', true);
    tc.verifyNotEmpty(h.peakMarkers);
    close(h.figure);
end

function testPlotScanBadScan(tc)
    bad = struct('counts', 1:10);   % no twoTheta
    tc.verifyError(@() xrdc.plot.plotScan(bad), 'xrdc:plot:badScan');
end

% ---------- publicationStyle ----------

function testPublicationStyleOnBlankAxes(tc)
    fig = figure();
    ax  = axes(fig);
    plot(ax, 1:10, (1:10).^2);
    xrdc.plot.publicationStyle(ax);   % should not error
    tc.verifyEqual(string(ax.YScale), "log");
    tc.verifyEqual(ax.FontSize, 18);
    close(fig);
end

function testPublicationStyleCustomSizes(tc)
    fig = figure();
    ax  = axes(fig);
    plot(ax, 1:10, 1:10);
    xrdc.plot.publicationStyle(ax, ...
        'TickFontSize', 14, 'LabelFontSize', 16, 'LogY', false);
    tc.verifyEqual(ax.FontSize, 14);
    tc.verifyEqual(string(ax.YScale), "linear");
    close(fig);
end

function testPublicationStyleRejectsNonAxes(tc)
    fig = figure();
    tc.verifyError(@() xrdc.plot.publicationStyle(fig), ...
        'xrdc:plot:badAxes');
    close(fig);
end

% ---------- plotStack ----------

function testPlotStackRendersAll(tc)
    scans = repmat(makeSyntheticScan(), 3, 1);
    scans(1).identifier = "a"; scans(2).identifier = "b"; scans(3).identifier = "c";
    h = xrdc.plot.plotStack(scans);
    tc.verifyLength(h.lines, 3);
    % Check the 3^(j) multiplier scaling: scan 2 max should be ~3x of scan 1 max
    y1 = h.lines(1).YData;
    y2 = h.lines(2).YData;
    tc.verifyEqual(max(y2) / max(y1), 3, 'RelTol', 1e-6);
    close(h.figure);
end

function testPlotStackCustomMultiplier(tc)
    scans = repmat(makeSyntheticScan(), 3, 1);
    h = xrdc.plot.plotStack(scans, 'Multiplier', 10);  % scalar base
    y1 = h.lines(1).YData;
    y3 = h.lines(3).YData;
    tc.verifyEqual(max(y3) / max(y1), 100, 'RelTol', 1e-6);
    close(h.figure);
end

function testPlotStackEmpty(tc)
    tc.verifyError(@() xrdc.plot.plotStack(repmat(makeSyntheticScan(), 0, 1)), ...
        'xrdc:plot:emptyStack');
end
