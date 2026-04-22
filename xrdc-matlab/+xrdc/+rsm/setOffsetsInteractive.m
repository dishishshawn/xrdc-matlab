function [deltaTheta, deltaOmega] = setOffsetsInteractive(scans, options)
%SETOFFSETSINTERACTIVE  Click-to-align RSM offsets against a known substrate peak.
%   [deltaTheta, deltaOmega] = xrdc.rsm.setOffsetsInteractive(scans)
%   [deltaTheta, deltaOmega] = xrdc.rsm.setOffsetsInteractive(scans, Name, Value, ...)
%
%   Implements the SetRSMOffsets workflow from xrdc1.pas:4020:
%     1. Plots the RSM with the current (or zero) offsets.
%     2. User clicks on the known substrate peak position.
%     3. The function computes domega for that slice.
%     4. A dialog asks for the theoretical 2θ and ω of the substrate peak.
%     5. Returns (ΔΘ, ΔΩ) in degrees.
%
%   The returned values can be passed directly to xrdc.rsm.toReciprocalSpace
%   or xrdc.plot.plotRsm as 'DeltaTheta' and 'DeltaOmega'.
%
%   Input
%     scans : (1×N) struct array from xrdc.rsm.loadAreaScan.
%
%   Name/Value
%     'Lambda'       (1,1) double  — wavelength in Å
%     'DeltaTheta'   (1,1) double  — initial ΔΘ for the preview plot (deg, default 0)
%     'DeltaOmega'   (1,1) double  — initial ΔΩ for the preview plot (deg, default 0)
%     'Flip'         (1,1) logical — negate k_par in preview (default false)
%
%   Output
%     deltaTheta  — 2θ correction in degrees
%     deltaOmega  — ω correction in degrees

    arguments
        scans                      (1,:) struct
        options.Lambda             (1,1) double  = NaN
        options.DeltaTheta         (1,1) double  = 0
        options.DeltaOmega         (1,1) double  = 0
        options.Flip               (1,1) logical = false
    end

    if isempty(scans)
        error('xrdc:rsm:emptyScans', 'scans array is empty.');
    end

    % Draw preview RSM so the user can see the current positions
    hFig = figure('Name', 'RSM — click the substrate peak', 'NumberTitle', 'off');
    xrdc.plot.plotRsm(scans, ...
        'Lambda',      options.Lambda, ...
        'DeltaTheta',  options.DeltaTheta, ...
        'DeltaOmega',  options.DeltaOmega, ...
        'Flip',        options.Flip, ...
        'TargetAxes',  gca);
    title(gca, 'Click on the substrate peak, then press Enter');

    % Capture the click in k-space coordinates
    [kParClick, kPerpClick] = ginput(1);
    title(gca, 'Offset dialog open — fill in the known peak position');

    % Invert k-space click back to approximate (2θ, ω) for the dialog default.
    % This is a rough estimate; the exact conversion is done from the user's input.
    lambda = options.Lambda;
    if isnan(lambda) && isfield(scans(1), 'lambda') && ~isnan(scans(1).lambda)
        lambda = scans(1).lambda;
    end

    if isnan(lambda)
        % Can't pre-fill; leave blanks
        approxTwoTheta = '';
        approxOmega    = '';
    else
        % Approximate: theta ~ asin(sqrt(kPar^2+kPerp^2)*lambda/2)
        k_total = sqrt(kParClick^2 + kPerpClick^2);
        if k_total > 0 && k_total * lambda / 2 <= 1
            theta_approx = asin(k_total * lambda / 2);   % radians
            tt_approx    = 2 * theta_approx * 180/pi;    % degrees
            % omega_approx: for symmetric reflection kPar≈0 → omega≈theta
            omega_approx = theta_approx * 180/pi + atan2(kParClick, kPerpClick)*180/pi;
            approxTwoTheta = sprintf('%.4f', tt_approx);
            approxOmega    = sprintf('%.4f', omega_approx);
        else
            approxTwoTheta = '';
            approxOmega    = '';
        end
    end

    % Dialog: user enters the known theoretical position of the substrate peak
    prompt  = {'Known substrate 2\theta (degrees):', ...
               'Known substrate \omega (degrees):'};
    dlgTitle = 'Set RSM Offsets — substrate peak position';
    defaults = {approxTwoTheta, approxOmega};
    answer   = inputdlg(prompt, dlgTitle, 1, defaults);

    if isempty(answer)
        close(hFig);
        error('xrdc:rsm:cancelled', 'setOffsetsInteractive cancelled by user.');
    end

    refTwoTheta = str2double(answer{1});
    refOmega    = str2double(answer{2});
    if isnan(refTwoTheta) || isnan(refOmega)
        close(hFig);
        error('xrdc:rsm:badInput', 'Could not parse the entered 2θ / ω values.');
    end

    % Find the slice with secondAxis closest to the clicked y-coordinate
    % (kPerp ~ ω in first approximation; but we match via secondAxis directly
    %  using the click's k_perp as a proxy for large-ω slices).
    % More robustly: find the slice whose k_perp centre is closest to kPerpClick.
    allSA = [scans.secondAxis];   % omega of each slice in degrees
    % 2θ_center of each slice (mean of first and last 2θ)
    allTTCtr = arrayfun(@(s) (s.twoTheta(1) + s.twoTheta(end)) / 2, scans);
    % domega = secondAxis - (2θ_center / 2) for each slice (degrees)
    domega_deg = allSA - allTTCtr / 2;

    % Identify which slice the user clicked on via the secondAxis closest
    % to the approximate omega inferred from the click
    if ~isnan(lambda) && ~isempty(approxOmega)
        [~, sliceIdx] = min(abs(allSA - str2double(approxOmega)));
    else
        [~, sliceIdx] = min(abs(allSA - refOmega));
    end

    % Compute ΔΘ and ΔΩ (xrdc1.pas:4020 logic)
    %   The measured peak is at (clickedTwoTheta, clickedOmega) in the current
    %   (un-offset) scan coordinates.  We want it to land at (refTwoTheta, refOmega).
    %
    %   From the ginput click in k-space, we recover the approximate 2θ and ω
    %   of the clicked point using the current (un-corrected) scan's domega:
    %       baseTwoTheta ≈ 2*arcsin(lambda*sqrt(kPar^2+kPerp^2)/2) * 180/pi
    %       baseOmega    ≈ baseTwoTheta/2 + domega_deg(sliceIdx)
    %
    %   deltaTheta = refTwoTheta - baseTwoTheta
    %   deltaOmega = refOmega    - baseOmega

    if ~isnan(lambda)
        k_total = sqrt(kParClick^2 + kPerpClick^2);
        if k_total > 0 && k_total * lambda / 2 <= 1
            theta_click  = asin(k_total * lambda / 2);           % rad
            baseTwoTheta = 2 * theta_click * 180/pi;             % deg
        else
            baseTwoTheta = refTwoTheta;  % fallback — zero correction
        end
    else
        baseTwoTheta = refTwoTheta;
    end
    baseOmega = baseTwoTheta / 2 + domega_deg(sliceIdx);

    deltaTheta = refTwoTheta - baseTwoTheta;
    deltaOmega = refOmega    - baseOmega;

    % Redraw with the computed offsets so the user can verify
    clf(hFig);
    xrdc.plot.plotRsm(scans, ...
        'Lambda',      lambda, ...
        'DeltaTheta',  deltaTheta, ...
        'DeltaOmega',  deltaOmega, ...
        'Flip',        options.Flip, ...
        'TargetAxes',  axes(hFig));
    title(gca, sprintf('Aligned RSM  \\DeltaΘ=%.4f°  \\DeltaΩ=%.4f°', ...
        deltaTheta, deltaOmega));
    drawnow;
end
