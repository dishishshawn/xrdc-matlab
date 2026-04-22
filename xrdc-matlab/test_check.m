addpath(pwd);
results = runtests;
for i = 1:numel(results)
    if ~results(i).Passed
        disp('---');
        disp(['Test: ' results(i).Name]);
        if isfield(results(i).Details, 'DiagnosticRecord')
            for k = 1:numel(results(i).Details.DiagnosticRecord)
                disp(results(i).Details.DiagnosticRecord(k).Message);
            end
        end
    end
end
