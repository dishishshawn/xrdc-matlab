%DEMOSTRAINRSM  RSM strain & composition on the real PtO2/TiO2 112 map.
%   Run from the xrdc-matlab root. Prints the analyzeStrainRSM report.
f = fullfile('..', 'data', ...
    'HP PtO2 on TiO2 001 112 RSM_C_HP PtO2 on TiO2 001 112 RSM_C.xrdml');
if ~isfile(f)
    error('demoStrainRSM:noData', 'Place the PtO2/TiO2 112 RSM in ../data.');
end
R = xrdc.rsm.analyzeStrainRSM(f, Substrate="TiO2", Film="PtO2", Reflection=[1 1 2]);
fprintf('Substrate %s: a=%.4f c=%.4f A (found=%d, off=%.2f%%)\n', ...
    R.substrate.name, R.substrate.aMeas, R.substrate.cMeas, ...
    R.substrate.found, R.substrate.predictOffPct);
if R.film.found
    fprintf('Film: a_par=%.4f a_perp=%.4f a0=%.4f A\n', R.aPar, R.aPerp, R.a0);
    fprintf('  strain par=%+.3f%% perp=%+.3f%%  relaxation=%.2f  pseudomorphic=%d\n', ...
        100*R.strainPar, 100*R.strainPerp, R.relaxation, R.pseudomorphic);
    if ~isnan(R.x), fprintf('  composition x=%.3f\n', R.x); end
else
    fprintf('Film peak not found.\n');
end
if ~isempty(R.flags), fprintf('Flags: %s\n', strjoin(cellstr(R.flags), ', ')); end
