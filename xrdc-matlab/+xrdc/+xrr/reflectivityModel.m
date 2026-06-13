function R = reflectivityModel(twoTheta, layers, lambda, options)
%REFLECTIVITYMODEL  Specular XRR via the Parratt recursion (Névot–Croce).
%   R = xrdc.xrr.reflectivityModel(twoTheta, layers, lambda)
%   R = xrdc.xrr.reflectivityModel(twoTheta, layers, lambda, Name=Value)
%
%   twoTheta : 2θ axis in degrees (the grazing angle is θ = 2θ/2).
%   layers   : struct array, top→bottom; the LAST element is the
%              semi-infinite substrate. Fields:
%                .material  (string)  material name/formula for opticalConstants
%                .density   (g/cm^3)
%                .thickness (nm)       finite for films; ignored for substrate
%                .roughness (nm)       roughness of this layer's TOP interface
%   lambda   : wavelength (Å). Cu Kα only (opticalConstants enforces this).
%
%   Name/Value
%     'Footprint'    (default true)  apply the low-angle illumination correction
%     'FootprintDeg' (default 0)     sample-fill angle in 2θ deg; 0 disables
%                                    the correction even when Footprint=true.
%
%   Returns R (Nx1) in [0,1], the specular reflectivity (no scale/background).
%   See also xrdc.xrr.opticalConstants, xrdc.xrr.fitReflectivity.

    arguments
        twoTheta            (:,1) double
        layers              (1,:) struct
        lambda              (1,1) double {mustBePositive}
        options.Footprint   (1,1) logical = true
        options.FootprintDeg(1,1) double {mustBeNonnegative} = 0
    end

    theta = twoTheta * (pi/180) / 2;          % grazing angle, rad
    k0    = 2*pi / lambda;                     % Å^-1
    L     = numel(layers);                     % media 1..L (last = substrate)

    % Complex refractive index per layer; ambient (vacuum) prepended.
    nAll = ones(L+1, 1) + 0i;                  % index 1 = ambient
    for j = 1:L
        [delta, beta] = xrdc.xrr.opticalConstants( ...
            string(layers(j).material), layers(j).density, lambda);
        nAll(j+1) = 1 - delta + 1i*beta;
    end

    % Vertical wavevector in each medium: kz = k0 * sqrt(n^2 - cos^2 θ).
    cos2 = cos(theta).^2;
    kz = zeros(numel(theta), L+1) + 0i;
    for m = 1:L+1
        kz(:,m) = k0 * sqrt(nAll(m)^2 - cos2);
    end

    % Parratt recursion from the substrate (X = 0) upward.
    % Media: 1=ambient, k+1 = layers(k). Interface between media k and k+1
    % has roughness layers(k).roughness; the layer BELOW it (media k+1 =
    % layers(k)) has thickness layers(k).thickness (substrate -> infinite).
    X = zeros(numel(theta), 1) + 0i;           % X at substrate top = 0
    for k = L:-1:1
        sigma = layers(k).roughness * 10;      % nm -> Å
        kzA = kz(:,k);  kzB = kz(:,k+1);
        r = (kzA - kzB) ./ (kzA + kzB) .* exp(-2 * kzA .* kzB * sigma^2);
        if k < L                                % media k+1 is a finite film
            d = layers(k).thickness * 10;       % nm -> Å
            phase = exp(2i * kzB * d);
        else                                    % substrate below: X already 0
            phase = 1;
        end
        X = (r + X .* phase) ./ (1 + r .* X .* phase);
    end

    R = abs(X).^2;

    % Footprint / illumination correction below the sample-fill angle.
    if options.Footprint && options.FootprintDeg > 0
        thetaFill = options.FootprintDeg * (pi/180) / 2;
        frac = min(1, sin(theta) ./ sin(thetaFill));
        frac(theta <= 0) = 0;
        R = R .* frac;
    end
end
