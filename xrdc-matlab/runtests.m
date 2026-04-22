function results = runtests()
%RUNTESTS  Run every xrdc unit-test suite and report.
%   Call from the xrdc-matlab root directory:
%       >> runtests
%
%   Returns the full matlab.unittest.TestResult array.

    here = fileparts(mfilename('fullpath'));
    addpath(here);                      % so +xrdc is visible
    suite = matlab.unittest.TestSuite.fromFolder(fullfile(here, 'tests'));
    runner = matlab.unittest.TestRunner.withTextOutput();
    results = runner.run(suite);
    disp(table(results));
end
