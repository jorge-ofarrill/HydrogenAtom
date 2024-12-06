function L = fast_laguerreL(n, k, x)
    % Computes the generalized Laguerre polynomial L_n^k(x)
    % n: Degree of the polynomial
    % k: Order of the polynomial
    % x: Evaluation points (vector or scalar)

    % Initialize for recurrence
    if n == 0
        L = ones(size(x)); % L_0^k(x) = 1
        return;
    elseif n == 1
        L = 1 + k - x; % L_1^k(x) = 1 + k - x
        return;
    end

    % Recurrence relation
    Lnm2 = ones(size(x)); % L_0^k(x)
    Lnm1 = 1 + k - x; % L_1^k(x)
    for m = 2:n
        L = ((2*m + k - 1 - x) .* Lnm1 - (m + k - 1) .* Lnm2) / m;
        Lnm2 = Lnm1; % Shift L_{n-2} to L_{n-1}
        Lnm1 = L; % Shift L_{n-1} to L_n
    end
end
