function [kPar, kPerp] = toReciprocalSpace(scan, options)
%TORERECIPROCALSPACE  Convert one θ-2θ slice to reciprocal-space coordinates.
%   [kPar, kPerp] = xrdc.rsm.toReciprocalSpace(scan)
%   [kPar, kPerp] = xrdc.rsm.toReciprocalSpace(scan, Name, Value, ...)
%
%   Implements the transform from ALGORITHM_SPEC §7.1 with the deliberate
%   θ-asymmetry from xrdc1.pas:3293-3296:
%
%     θ_raw  = 2θ_point / 2            (uncorrected; used to build ω)
%     ω      = secondAxis − (2θ_ctr/2) + θ_raw + ΔΩ   (radians)
%     θ      = (2θ_point + ΔΘ) / 2    (corrected; used in k formulas)
%     k_perp = (2/λ) · sin(θ) · cos(ω − θ)
%     k_par  = ±(2/λ) · sin(θ) · sin(ω − θ)
%
%   where 2θ_ctr = mean([2θ(1), 2θ(end)]) normalises the ω baseline so
%   that all slices in an area scan share a consistent k-space origin.
%
%   For non-twoThetaOmega scans the simplified branch is used:
%     ω = (secondAxis + ΔΩ) · π/180
%     θ = (2θ_point  + ΔΘ) · π/360
%
%   WARNING: do not "fix" the θ asymmetry. When ΔΘ ≠ 0 the two uses of θ
%   intentionally differ. Undoing this shifts RSM peaks off-target.
%
%   Input
%     scan        : struct with .twoTheta (deg), .secondAxis (deg),
%                   .scanType, and .lambda (Å, or NaN).
%
%   Name/Value
%     'Lambda'       (1,1) double  — wavelength in Å; overrides scan.lambda.
%                                    Required if scan.lambda is NaN.
%     'DeltaTheta'   (1,1) double  — 2θ offset in degrees  (ΔΘ, default 0)
%     'DeltaOmega'   (1,1) double  — ω offset in degrees   (ΔΩ, default 0)
%     'Flip'         (1,1) logical — negate k_par (default false)
%
%   Output
%     kPar   (Nx1 double) — in-plane reciprocal coordinate, Å⁻¹
%     kPerp  (Nx1 double) — out-of-plane reciprocal coordinate, Å⁻¹

    arguments
        scan                     (1,1) struct
        options.Lambda           (1,1) double  = NaN
        options.DeltaTheta       (1,1) double  = 0
        options.DeltaOmega       (1,1) double  = 0
        options.Flip             (1,1) logical = false
    end

    % Resolve wavelength
    lambda = options.Lambda;
    if isnan(lambda)
        if isfield(scan, 'lambda') && ~isnan(scan.lambda)
            lambda = scan.lambda;
        else
            error('xrdc:rsm:noLambda', ...
                'Wavelength not found in scan.lambda; supply ''Lambda'' option.');
        end
    end
    if lambda <= 0
        error('xrdc:rsm:badLambda', 'Lambda must be positive (in Å).');
    end

    if ~isfield(scan, 'twoTheta') || isempty(scan.twoTheta)
        error('xrdc:rsm:badScan', 'scan.twoTheta is missing or empty.');
    end
    if ~isfield(scan, 'secondAxis') || isnan(scan.secondAxis)
        error('xrdc:rsm:noSecondAxis', ...
            'scan.secondAxis is NaN; required for RSM transform.');
    end

    tt_deg  = double(scan.twoTheta(:));           % 2θ in degrees
    dTheta  = options.DeltaTheta;                 % degrees
    dOmega  = options.DeltaOmega;                 % degrees
    sa_deg  = double(scan.secondAxis);            % secondAxis in degrees

    scanType = "twoThetaOmega";
    if isfield(scan, 'scanType')
        scanType = scan.scanType;
    end

    if scanType == "twoThetaOmega"
        % Full asymmetric branch (xrdc1.pas:3223-3296)
        tt_center_deg = (tt_deg(1) + tt_deg(end)) / 2;  % 2θ_ctr in degrees

        % Convert all angles to radians for computation
        tt_rad        = tt_deg        * (pi/180);
        tt_center_rad = tt_center_deg * (pi/180);
        sa_rad        = sa_deg        * (pi/180);
        dTheta_rad    = dTheta        * (pi/180);
        dOmega_rad    = dOmega        * (pi/180);

        % θ_raw: uncorrected half-angle — used ONLY to build ω
        theta_raw = tt_rad / 2;

        % ω: uses uncorrected θ_raw (the asymmetry — do not change)
        omega = sa_rad - (tt_center_rad / 2) + theta_raw + dOmega_rad;

        % θ: corrected half-angle — used in k formulas
        theta = (tt_rad + dTheta_rad) / 2;
    else
        % Simplified branch for non-twoThetaOmega scans
        dTheta_rad = dTheta * (pi/180);
        dOmega_rad = dOmega * (pi/180);
        omega = (sa_deg  + dOmega) * (pi/180);
        theta = (tt_deg  + dTheta) * (pi/360);
    end

    scale = 2 / lambda;
    sinT  = sin(theta);
    diff  = omega - theta;

    kPerp =  scale .* sinT .* cos(diff);
    kPar  =  scale .* sinT .* sin(diff);

    if options.Flip
        kPar = -kPar;
    end
end
