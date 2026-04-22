addpath(pwd);
x = (0:0.01:10).';
y = x.^2;
[~, slope2] = xrdc.signal.derivatives(x, y, 15, 3);
mid = 50:numel(x)-50;

disp('Expected (should be all 2):');
disp(slope2(mid(1:10)));
disp(' ');
disp('Min:');
disp(min(slope2(mid)));
disp('Max:');
disp(max(slope2(mid)));
disp('Mean:');
disp(mean(slope2(mid)));
