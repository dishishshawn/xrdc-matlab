%DEMOXRRFIT  Slab-model XRR fitting: thickness, density, roughness.
%   Loads an XRR scan, fits a single film on a substrate, prints the
%   parameters with uncertainties, and overlays the model on the data.
%
%   Usage: edit `path`, `film`, `sub` for your scan, then run.

path = fullfile('..','..','data', ...
    'TR_S25_SRO_STO(100)_700c_100mT_10500sh_5hz_XRR_05062026.txt');
film = "SrRuO3"; sub = "SrTiO3";

scan = xrdc.io.readScan(path);
res  = xrdc.xrr.fitReflectivity(scan, Film=film, Substrate=sub);

fprintf('Thickness : %.2f ± %.2f nm\n', res.thicknessNm, res.thicknessSeNm);
fprintf('Density   : %.2f g/cm^3 (%.0f%% of bulk)\n', ...
    res.densityGcc, 100*res.densityFraction);
fprintf('Roughness : %.2f nm film / %.2f nm interface\n', ...
    res.filmRoughnessNm, res.substrateRoughnessNm);
fprintf('R^2 = %.4f  (%s)\n', res.rSquared, res.method);

figure;
semilogy(scan.twoTheta, max(scan.counts,1), '.', 'DisplayName','data'); hold on;
semilogy(scan.twoTheta, max(res.modelCurve,1), '-', 'LineWidth',1.5, ...
    'DisplayName','Parratt fit');
xlabel('2\theta (°)'); ylabel('Counts'); legend; title('XRR slab-model fit');
