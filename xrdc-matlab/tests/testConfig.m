function tests = testConfig
%TESTCONFIG  Unit tests for the +xrdc.+config lab-preset package.
    tests = functiontests(localfunctions);
end

function tmp = tempJson()
%TEMPJSON  A unique, non-existent JSON path under a temp folder.
    tmp = string(fullfile(tempname, 'settings.json'));
end

% ---------------------------------------------------------------------
function testDefaultsShape(tc)
    cfg = xrdc.config.defaults();
    tc.verifyClass(cfg, 'struct');
    tc.verifyEqual(sort(fieldnames(cfg)), sort({ ...
        'wavelength'; 'applyWavelength'; 'substrate'; 'rcShape'; ...
        'useFilenameRules'; 'defaultScanType'}));
    % Defaults reproduce the historical Paik-lab Rigaku behaviour.
    tc.verifyEqual(cfg.wavelength, 1.5406, 'AbsTol', 1e-9);
    tc.verifyTrue(cfg.applyWavelength);
    tc.verifyEqual(cfg.substrate, "SrTiO3");
    tc.verifyEqual(cfg.rcShape, "gauss");
    tc.verifyTrue(cfg.useFilenameRules);
    tc.verifyEqual(cfg.defaultScanType, "twoThetaOmega");
end

function testDefaultsTypes(tc)
    cfg = xrdc.config.defaults();
    tc.verifyClass(cfg.wavelength, 'double');
    tc.verifyClass(cfg.applyWavelength, 'logical');
    tc.verifyClass(cfg.substrate, 'string');
    tc.verifyClass(cfg.rcShape, 'string');
    tc.verifyClass(cfg.useFilenameRules, 'logical');
    tc.verifyClass(cfg.defaultScanType, 'string');
end

function testLoadMissingReturnsDefaults(tc)
    cfg = xrdc.config.load(tempJson());
    tc.verifyEqual(cfg, xrdc.config.defaults());
end

function testSaveLoadRoundTrip(tc)
    p = tempJson();
    c = onCleanup(@() cleanupPath(p));

    cfg = xrdc.config.defaults();
    cfg.wavelength       = 0.7093;          % Mo Kα1
    cfg.applyWavelength  = false;
    cfg.substrate        = "LaAlO3";
    cfg.rcShape          = "lorentz";
    cfg.useFilenameRules = false;
    cfg.defaultScanType  = "omega";

    xrdc.config.save(cfg, p);
    tc.verifyTrue(isfile(p));

    back = xrdc.config.load(p);
    tc.verifyEqual(back.wavelength, 0.7093, 'AbsTol', 1e-9);
    tc.verifyFalse(back.applyWavelength);
    tc.verifyEqual(back.substrate, "LaAlO3");
    tc.verifyEqual(back.rcShape, "lorentz");
    tc.verifyFalse(back.useFilenameRules);
    tc.verifyEqual(back.defaultScanType, "omega");
end

function testSaveCreatesFolder(tc)
    % save() must create a missing parent folder on demand.
    p = string(fullfile(tempname, 'nested', 'deep', 'settings.json'));
    c = onCleanup(@() cleanupPath(p));
    xrdc.config.save(xrdc.config.defaults(), p);
    tc.verifyTrue(isfile(p));
end

function testLoadMergesPartialFile(tc)
    % A file with only some fields keeps defaults for the rest.
    p = tempJson();
    c = onCleanup(@() cleanupPath(p));
    mkdir(fileparts(p));
    writeText(p, '{"substrate":"TiO2","wavelength":1.79}');

    cfg = xrdc.config.load(p);
    tc.verifyEqual(cfg.substrate, "TiO2");
    tc.verifyEqual(cfg.wavelength, 1.79, 'AbsTol', 1e-9);
    % Untouched fields fall back to defaults.
    tc.verifyEqual(cfg.rcShape, "gauss");
    tc.verifyTrue(cfg.useFilenameRules);
end

function testLoadIgnoresUnknownFields(tc)
    p = tempJson();
    c = onCleanup(@() cleanupPath(p));
    mkdir(fileparts(p));
    writeText(p, '{"substrate":"TiO2","bogusField":42}');

    cfg = xrdc.config.load(p);
    tc.verifyEqual(cfg.substrate, "TiO2");
    tc.verifyFalse(isfield(cfg, 'bogusField'));
end

function testSaveDropsExtraFields(tc)
    % An over-stuffed struct only writes the canonical fields.
    p = tempJson();
    c = onCleanup(@() cleanupPath(p));
    cfg = xrdc.config.defaults();
    cfg.extra = "should not persist";
    xrdc.config.save(cfg, p);

    back = xrdc.config.load(p);
    tc.verifyFalse(isfield(back, 'extra'));
    tc.verifyEqual(sort(fieldnames(back)), sort(fieldnames(xrdc.config.defaults())));
end

function testLoadBadJsonWarnsAndDefaults(tc)
    p = tempJson();
    c = onCleanup(@() cleanupPath(p));
    mkdir(fileparts(p));
    writeText(p, '{ this is not valid json ');

    cfg = tc.verifyWarning(@() xrdc.config.load(p), 'xrdc:config:badFile');
    tc.verifyEqual(cfg, xrdc.config.defaults());
end

function testCoerceTypesFromJson(tc)
    % Logicals stored as JSON true/false and numbers as JSON numbers come
    % back with the documented classes.
    p = tempJson();
    c = onCleanup(@() cleanupPath(p));
    mkdir(fileparts(p));
    writeText(p, '{"applyWavelength":false,"useFilenameRules":true,"wavelength":2.29}');

    cfg = xrdc.config.load(p);
    tc.verifyClass(cfg.applyWavelength, 'logical');
    tc.verifyFalse(cfg.applyWavelength);
    tc.verifyTrue(cfg.useFilenameRules);
    tc.verifyClass(cfg.wavelength, 'double');
    tc.verifyEqual(cfg.wavelength, 2.29, 'AbsTol', 1e-9);
end

function testConfigPathUnderPrefdir(tc)
    p = xrdc.config.configPath();
    tc.verifyClass(p, 'string');
    tc.verifyTrue(startsWith(p, string(prefdir)));
    tc.verifyTrue(endsWith(p, "settings.json"));
end

% ---------------------------------------------------------------------
function writeText(path, txt)
    fid = fopen(char(path), 'w');
    assert(fid >= 0, 'could not open %s', path);
    cleanup = onCleanup(@() fclose(fid));
    fwrite(fid, txt, 'char');
end

function cleanupPath(path)
    folder = fileparts(path);
    if isfolder(folder)
        rmdir(char(folder), 's');
    elseif isfile(path)
        delete(char(path));
    end
end
