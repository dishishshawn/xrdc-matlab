% Quick test of derivatives fallback
addpath(pwd);
try
    x = (0:0.01:10).';
    y = 3*x;
    disp('Testing derivatives with fallback...');
    [slope, slope2] = xrdc.signal.derivatives(x, y, 11, 3);
    disp('Success! Slope sample:');
    disp(slope(20:25));
catch e
    disp('ERROR:');
    disp(e.message);
    disp(e.stack);
end
