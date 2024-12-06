function R_nl = radial_wavefunction(n, l, r, a0)
    % Radial wavefunction for hydrogen atom
    % n: principal quantum number
    % l: azimuthal quantum number
    % r: radial distance (same as rho in z=0 plane)
    % a0: Bohr radius

    % Normalization factor
    normalization = sqrt((2 / (n * a0))^3 * factorial(n - l - 1) / ...
                         (2 * n * factorial(n + l)^3));

    % Dimensionless variable rho = 2r / (n * a0)
    rho = 2 * r / (n * a0);

    % Associated Laguerre polynomial
    L = fast_laguerreL(n - l - 1, 2 * l + 1, rho);

    % Radial wavefunction
    R_nl = normalization * rho.^l .* exp(-rho / 2) .* L;
end
