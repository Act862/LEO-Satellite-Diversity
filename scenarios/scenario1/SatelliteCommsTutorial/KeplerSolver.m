function E = KeplerSolver(M, e)
    % Solve Kepler's Equation: M = E - e*sin(E)
    % Use Newton-Raphson Method
    E = M; % Initial guess
    tol = 1e-8;
    max_iter = 1000;
    iter = 0;

    while iter < max_iter
        f = E - e * sin(E) - M;
        f_prime = 1 - e * cos(E);
        E_new = E - f / f_prime;

        if abs(E_new - E) < tol
            break;
        end

        E = E_new;
        iter = iter + 1;
    end
end
