function Y_lm = spherical_harmonic(l, m, theta, phi)
    % Spherical harmonic for hydrogen atom
    % l: azimuthal quantum number
    % m: magnetic quantum number
    % theta: polar angle (radians)
    % phi: azimuthal angle (radians)

    % Normalization factor
    normalization = sqrt((2 * l + 1) / (4 * pi) * factorial(l - abs(m)) / factorial(l + abs(m)));

    % Compute associated Legendre polynomial using built-in function
    % theta = gpuArray(theta);
    P_lm_array = legendre(l, cos(theta), 'norm'); % Normalized associated Legendre
    P_lm = squeeze(P_lm_array(abs(m) + 1, :,:,:)); % Extract the correct m

    % Adjust for m < 0 (factor of (-1)^m for negative m values)
    if m < 0
        P_lm = P_lm * (-1)^m;
    end

    % Complex exponential for azimuthal dependence
    exp_imphi = exp(1i * m * phi);

    % Spherical harmonic
    Y_lm = normalization .* P_lm .* exp_imphi;
    % keyboard
end
