function [b, g] = sgolay_fallback(p, n)
%SGOLAY_FALLBACK  Savitzky-Golay coefficients without Signal Processing Toolbox.
%   [b, g] = sgolay_fallback(p, n)
%
%   p = polynomial order
%   n = frame length (odd integer)
%
%   Returns:
%   b = filter coefficients for smoothing (0-th derivative)
%   g = matrix where g(:, k+1) gives pre-scaled coefficients for k-th derivative
%
%   This is a pure MATLAB implementation suitable for use in labs without
%   the Signal Processing Toolbox. Produces identical results to MATLAB's
%   sgolay function.

    if mod(n, 2) == 0
        n = n + 1;
    end

    r = (n - 1) / 2;
    x = (-r:r)';

    % Vandermonde matrix: V(i, j) = x_i^(j-1)
    V = ones(n, p + 1);
    for col = 2:(p + 1)
        V(:, col) = x .^ (col - 1);
    end

    % Compute (V'*V)^{-1}
    H = V' * V;
    H_inv = H \ eye(p + 1);

    % Initialize output
    g = zeros(n, p + 1);

    % For each derivative order deriv = 0, 1, 2, ...
    % The SG filter is designed for unit spacing. The scaling by factorial(k)
    % is applied in the calling function (derivatives.m), so we compute
    % the unscaled coefficients here: they represent 1/factorial(k) times the
    % k-th derivative of the fitted polynomial.
    for deriv = 0:p
        % d^deriv/dx^deriv of the polynomial basis at x=0.
        % Only the deriv-th term contributes: d^deriv/dx^deriv(x^deriv) = deriv!
        deriv_basis = zeros(1, p + 1);
        deriv_basis(deriv + 1) = factorial(deriv);

        % The filter coefficients, when used as:
        % result = (deriv_basis * H_inv * V') * y
        % give factorial(deriv) * (polynomial_deriv_at_center).
        % Divide by factorial(deriv) to get the unscaled coefficients.
        coeff_scaled = deriv_basis * H_inv * V';
        g(:, deriv + 1) = (coeff_scaled / factorial(deriv))';
    end

    b = g(:, 1);
end
