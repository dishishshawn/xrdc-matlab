function tests = testIo
%TESTIO  Unit tests for the +xrdc.+io package.
%   Limited until real sample files land; these tests exercise the text
%   reader and the dispatcher logic against temporary files.
    tests = functiontests(localfunctions);
end

function testReadTextScanBasic(tc)
    tmp = [tempname '.txt'];
    cleanup = onCleanup(@() delete(tmp));
    twoTheta = (20:0.01:30).';
    counts   = round(1e4 * exp(-((twoTheta - 25).^2)/0.5));
    writematrix([twoTheta, counts], tmp, 'Delimiter', 'tab', 'FileType', 'text');
    scan = xrdc.io.readTextScan(tmp);
    tc.verifyEqual(scan.twoTheta, twoTheta, 'AbsTol', 1e-6);
    tc.verifyEqual(scan.counts,   counts,   'AbsTol', 1e-6);
    tc.verifyEqual(scan.sourceFormat, "text");
end

function testReadTextScanCommaDecimal(tc)
    % German locale — commas as decimal separators
    tmp = [tempname '.txt'];
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w');
    fprintf(fid, '20,00\t1000\n');
    fprintf(fid, '20,01\t1050\n');
    fprintf(fid, '20,02\t1100\n');
    fclose(fid);
    scan = xrdc.io.readTextScan(tmp);
    tc.verifyEqual(scan.twoTheta, [20.00; 20.01; 20.02], 'AbsTol', 1e-6);
    tc.verifyEqual(scan.counts,   [1000; 1050; 1100],    'AbsTol', 1e-6);
end

function testReadTextScanSkipsComments(tc)
    tmp = [tempname '.txt'];
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w');
    fprintf(fid, '# Sample scan header\n');
    fprintf(fid, '%% MATLAB-style comment\n');
    fprintf(fid, 'two_theta  counts\n');   % non-numeric header
    fprintf(fid, '20.0  100\n');
    fprintf(fid, '20.1  200\n');
    fclose(fid);
    scan = xrdc.io.readTextScan(tmp);
    tc.verifyEqual(scan.twoTheta, [20.0; 20.1], 'AbsTol', 1e-10);
    tc.verifyEqual(scan.counts,   [100; 200]);
end

function testDispatcherTextFallback(tc)
    tmp = [tempname '.dat'];
    cleanup = onCleanup(@() delete(tmp));
    writematrix([10 100; 11 200; 12 300], tmp, 'Delimiter','tab', 'FileType','text');
    scan = xrdc.io.readScan(tmp);
    tc.verifyEqual(scan.sourceFormat, "text");
    tc.verifyEqual(scan.counts, [100; 200; 300]);
end

function testDispatcherRigakuBlocked(tc)
    tmp = [tempname '.raw'];
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w');
    fprintf(fid, 'RAW1.01 fake content\n');
    fclose(fid);
    tc.verifyError(@() xrdc.io.readScan(tmp), 'xrdc:io:notImplemented');
end

function testEmptyScanShape(tc)
    s = xrdc.io.emptyScan();
    expectedFields = {'twoTheta','counts','scanType','secondAxis', ...
                      'secondAxisName','identifier','lambda','metadata', ...
                      'sourcePath','sourceFormat'};
    tc.verifyEqual(sort(fieldnames(s)), sort(expectedFields.'));
end

% =====================================================================
% Rigaku SmartLab .txt export
% =====================================================================

function testReadRigakuTxtHeadered(tc)
    % Synthetic file mimicking the Paik lab's Rigaku SmartLab export.
    tmp = [tempname 'TR_synthetic_2theta_omega.txt'];
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w', 'n', 'UTF-8');
    fwrite(fid, char([239 187 191]));   % UTF-8 BOM
    fprintf(fid, 'TR_synthetic_2theta_omega\t\n');
    fprintf(fid, '2θ, °\tIntensity, cps\n');
    for v = 20:0.01:20.1
        fprintf(fid, '%.6f\t%.4f\n', v, 1000 + 10000*exp(-((v-20.05)^2)/0.001));
    end
    fclose(fid);

    scan = xrdc.io.readRigakuTxt(tmp);
    tc.verifyEqual(scan.sourceFormat, "rigakuTxt");
    tc.verifyEqual(scan.scanType,     "twoThetaOmega");
    tc.verifyEqual(numel(scan.twoTheta), 11);
    tc.verifyEqual(scan.twoTheta(1),  20.00, 'AbsTol', 1e-6);
    tc.verifyEqual(scan.twoTheta(end),20.10, 'AbsTol', 1e-6);
    tc.verifyTrue(scan.metadata.headered);
end

function testReadRigakuTxtHeaderless(tc)
    % Some Rigaku XRR exports arrive without the 2-line header.
    tmp = [tempname 'TR_synthetic_XRR.txt'];
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w');
    for v = 0:0.01:0.05
        fprintf(fid, '%.6f\t%.4f\n', v, 2e6);
    end
    fclose(fid);

    scan = xrdc.io.readRigakuTxt(tmp);
    tc.verifyEqual(scan.sourceFormat, "rigakuTxt");
    tc.verifyEqual(scan.scanType,     "twoThetaOmega");   % XRR = coupled
    tc.verifyEqual(numel(scan.twoTheta), 6);
    tc.verifyFalse(scan.metadata.headered);
end

function testReadRigakuTxtRockingCurveInference(tc)
    % Filename contains "RC" → scanType should be "omega".
    tmp = [tempname 'TR_sample_film_RC_date.txt'];
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w');
    fprintf(fid, 'TR_sample_film_RC_date\n');
    fprintf(fid, '2θ, °\tIntensity, cps\n');
    for v = -0.5:0.05:0.5
        fprintf(fid, '%.4f\t%.1f\n', v, 100 + 1000*exp(-v^2/0.01));
    end
    fclose(fid);

    scan = xrdc.io.readRigakuTxt(tmp);
    tc.verifyEqual(scan.scanType, "omega");
end

function testReadRigakuTxtPhiInference(tc)
    tmp = [tempname 'TR_sample_101_phi_scan.txt'];
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w');
    fprintf(fid, 'TR_sample_101_phi_scan\n');
    fprintf(fid, '2θ, °\tIntensity, cps\n');
    for v = 0:1:359
        fprintf(fid, '%d\t%.1f\n', v, 50);
    end
    fclose(fid);

    scan = xrdc.io.readRigakuTxt(tmp);
    tc.verifyEqual(scan.scanType, "phi");
end

function testDispatcherRoutesRigakuTxt(tc)
    % readScan must pick readRigakuTxt via the "INTENSITY, CPS" marker.
    tmp = [tempname 'TR_dispatcher_test.txt'];
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w');
    fprintf(fid, 'TR_dispatcher_test\n');
    fprintf(fid, '2θ, °\tIntensity, cps\n');
    fprintf(fid, '20.0\t1000\n20.1\t1100\n');
    fclose(fid);
    scan = xrdc.io.readScan(tmp);
    tc.verifyEqual(scan.sourceFormat, "rigakuTxt");
end

function testDispatcherTrPrefixRoutesToRigaku(tc)
    % Even without the "INTENSITY, CPS" marker, the TR_ filename prefix
    % triggers the Rigaku parser (headerless XRR variant).
    tmp = fullfile(tempdir, 'TR_headerless_XRR.txt');
    cleanup = onCleanup(@() delete(tmp));
    fid = fopen(tmp, 'w');
    fprintf(fid, '0.0\t2000000\n0.01\t2100000\n0.02\t2200000\n');
    fclose(fid);
    scan = xrdc.io.readScan(tmp);
    tc.verifyEqual(scan.sourceFormat, "rigakuTxt");
end

function testReadRigakuTxtRealFiles(tc)
    % Integration: read every .txt file in the Paik lab Rigaku drop.
    dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), ...
        '..', '..', 'rexdrctomatlabport_rigakudatasets');
    assumeTrue(tc, isfolder(dataDir), ...
        'Rigaku data folder not present; skipping real-file test.');

    files = dir(fullfile(dataDir, '*.txt'));
    assumeFalse(tc, isempty(files), 'No .txt files found in Rigaku data folder.');

    for i = 1:numel(files)
        p = fullfile(files(i).folder, files(i).name);
        scan = xrdc.io.readScan(p);
        tc.verifyEqual(scan.sourceFormat, "rigakuTxt", ...
            sprintf('Dispatcher misrouted %s', files(i).name));
        tc.verifyGreaterThan(numel(scan.twoTheta), 10, ...
            sprintf('Too few points in %s', files(i).name));
        tc.verifyTrue(all(isfinite(scan.twoTheta)), ...
            sprintf('Non-finite 2θ in %s', files(i).name));
        tc.verifyTrue(all(isfinite(scan.counts)), ...
            sprintf('Non-finite counts in %s', files(i).name));
    end
end
