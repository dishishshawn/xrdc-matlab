addpath(pwd);

% Test sgolay_fallback vs expected values for 2nd deriv
frameSize = 3;
polyOrder = 2;

[b, g_fallback] = xrdc.signal.sgolay_fallback(polyOrder, frameSize);

disp('sgolay_fallback for frameSize=3, polyOrder=2:');
disp('0-th deriv (smoothing):');
disp(g_fallback(:, 1)');
disp('1st deriv:');
disp(g_fallback(:, 2)');
disp('2nd deriv:');
disp(g_fallback(:, 3)');

% Test on y = [1, 0, 1] (which is x^2 at x=[-1,0,1])
y = [1; 0; 1];
approx_2nd_deriv = 2 * (g_fallback(:, 3)' * y) / 1^2;
disp(' ');
disp('For y = [1, 0, 1] (x^2 at unit spacing):');
disp(['Computed 2nd deriv: ' num2str(approx_2nd_deriv)]);
disp('Expected: 2');

% Compare with real sgolay if available
if ~isempty(which('sgolay'))
    [b_real, g_real] = sgolay(polyOrder, frameSize);
    disp(' ');
    disp('Real sgolay 2nd deriv coeffs:');
    disp(g_real(:, 3)');
    approx_2nd_deriv_real = 2 * (g_real(:, 3)' * y) / 1^2;
    disp(['Computed 2nd deriv with real sgolay: ' num2str(approx_2nd_deriv_real)]);
end
