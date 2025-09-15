%% Doppler Fading Function (Jakes' model approximation)
function H = doppler_fading(numSymbols, fd, M, K)
    if nargin < 4
        K = 0; % Pure Rayleigh if K not provided
    end
    t = (0:numSymbols-1);
    H = zeros(M, numSymbols);
    for m = 1:M
        phi = 2 * pi * rand(1, 8);
        f = fd * cos(2 * pi * (0:7) / 8);
        Jakes_wave = sqrt(2 / 8) * sum(cos(2 * pi * f.' * t + phi.'), 1);
        if K > 0
            % Rician fading with LOS component
            LOS = sqrt(K / (K + 1));
            NLOS = sqrt(1 / (2 * (K + 1)));
            H(m, :) = LOS + NLOS * Jakes_wave .* exp(1j * 2 * pi * rand(1, numSymbols));
        else
            % Pure Rayleigh fading
            H(m, :) = Jakes_wave .* exp(1j * 2 * pi * rand(1, numSymbols));
        end
    end
end