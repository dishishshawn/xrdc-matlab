%% Demo: predict SrTiO3 substrate peaks and compare to the tabulated values.
% Run from the xrdc-matlab root so the +xrdc package is on the path.

addpath(fileparts(fileparts(mfilename('fullpath'))));

%% 1. Wavelength from Cu Kα1 (default)
energyEv = 8049.19;
lambda = xrdc.lattice.energyToLambda(energyEv);
fprintf('Cu Kα1 wavelength: %.4f Å (from E=%.2f eV)\n', lambda, energyEv);

%% 2. Simulate SrTiO3 pattern for hkl ∈ [0, 4]
lat = struct('system', 'cubic', 'a', 3.905);
T = xrdc.lattice.simulatePattern(lat, [0 4 0 4 0 4], lambda, ...
    'TwoThetaRange', [20, 120]);

%% 3. Compare to tabulated peaks from substrates.json
substratesPath = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
    '+xrdc', '+data', 'substrates.json');
sub = jsondecode(fileread(substratesPath));
tabulated = sub.substrates.SrTiO3_x00_.peaks;   % JSON key got sanitised

fprintf('\nhkl    predicted 2θ   tabulated 2θ   Δ\n');
for hkl = ["100", "200", "300", "400"]
    idx = find(T.label == "(" + extractBefore(hkl,2) + " " + extractBetween(hkl,2,2) + " " + extractAfter(hkl,2) + ")", 1);
    if isempty(idx), continue, end
    tab = tabulated.("x" + hkl);   % jsondecode prefixes all-digit keys with 'x'
    fprintf('%s    %8.3f°     %8.3f°   %+6.3f°\n', ...
        hkl, T.twoTheta(idx), tab, T.twoTheta(idx) - tab);
end

%% 4. Round-trip the Nelson-Riley workflow on a synthetic Si calibrant
% Use the tabulated Si peaks from Substrates.def with their assumed m (order)
lambda = 1.5406;
siPeaks = [28.442, 47.302, 56.121, 69.130, 76.377, 88.026, 94.948, ...
           106.715, 114.087, 127.541, 136.890];
hklSi   = {'111','220','311','400','331','422','511','440','531','620','533'};
% compute per-peak "a" from Bragg's law and the known order family
% (assume cubic Si, a≈5.43088 Å)
aExpected = 5.43088;
d = xrdc.lattice.twoThetaToD(siPeaks, lambda);

% m_i = sqrt(h²+k²+l²) for cubic
sumSq = zeros(size(hklSi));
for i = 1:numel(hklSi)
    hkl = hklSi{i};
    sumSq(i) = sum((hkl - '0').^2);   % '1'->1, '0'->0 etc.
end
aPerPeak = d .* sqrt(sumSq);

nr = xrdc.lattice.nelsonRiley(siPeaks.', aPerPeak.');
fprintf('\nNelson-Riley Si calibration:\n');
fprintf('  a₀ = %.5f ± %.5f Å   (R² = %.4f)\n', nr.a0, nr.a0SE, nr.rSquared);
fprintf('  expected 5.43088 Å\n');
