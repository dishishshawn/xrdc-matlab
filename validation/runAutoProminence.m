%RUNAUTOPROMINENCE  Auto-prominence sweep over the theta-2theta data dump.
%   Prints, per scan: peak count, then 2-theta / counts / FWHM per peak.
%   Eyeball criteria (design spec 2026-06-12, conservative threshold):
%   STO scans show the 3 substrate peaks (22.75/46.47/72.57 +/- 0.05,
%   except S10 which is sample-shifted ~0.1 deg) plus any CONFIDENT film
%   peaks; weak ~3x-background film bumps are correctly omitted; no two
%   peaks closer than ~0.15 deg; no baseline-noise picks.
%   Run from the repo root with xrdc-matlab on the path:
%     cd xrdc-matlab && matlab -batch "addpath(pwd); cd ../validation; runAutoProminence"

names = { ...
    'TR_S04_PTO_STO(100)_750c_200mT_1000sh_3hz_2theta omega_04072026.txt', ...
    'TR_S05_PTO_STO(100)_600c_200mT_1000sh_2hz_2theta omega_04092026.txt', ...
    'TR_S06_PTO_LAO(100)_600c_200mT_1000sh_2hz_2 theta omega_04082026.txt', ...
    'TR_S07_PTO_STO(100)_600c_variable pressure_2000sh_3hz_2theta omega_04102026.txt', ...
    'TR_S08_PTO_STO(100)_550c_200mT_5000sh_3hz_2theta omega_04132026.txt', ...
    'TR_S10_PTO_STO(100)_500c_150mT_20000sh_5hz_2theta omega_04162026.txt', ...
    'TR_S11_PTO_STO(100)_580c_150mT_20000sh_5hz_2theta omega_04162026.txt'};

dataDir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
for k = 1:numel(names)
    p = fullfile(dataDir, names{k});
    if ~isfile(p)
        fprintf('\nSKIP (missing): %s\n', names{k});
        continue
    end
    scan = xrdc.io.readScan(p);
    pk = xrdc.peaks.findPeaks(scan, 'MinSeparation', 0.2);
    fprintf('\n%s\n  -> %d peak(s)\n', names{k}, numel(pk));
    for i = 1:numel(pk)
        fprintf('   %8.3f deg   %10.0f cts   fwhm %.3f deg\n', ...
            pk(i).twoTheta, pk(i).counts, pk(i).fwhm);
    end
end
