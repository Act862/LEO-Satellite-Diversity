% function h = loglogistic_fading(M, N, shape, scale)
% % Generates log-logistic fading amplitudes for M antennas and N samples
% % shape = shape parameter (alpha > 0)
% % scale = scale parameter (beta > 0)
% % Formula for CDF: F(x) = 1 / (1 + (x/scale)^(-shape)) for x > 0
% % Generate by inverse transform sampling:
% U = rand(M, N); % uniform(0,1)
% h = scale * (U ./ (1 - U)).^(1/shape);
% % For complex channel, multiply by uniform random phase
% phase = exp(1j * 2 * pi * rand(M, N));
% h = h .* phase;
% end
%% Log-logistic Fading Function
function H = loglogistic_fading(M, numSymbols, alpha, beta)
    U = rand(M, numSymbols);
    magnitude = beta * (U ./ (1 - U)).^(1 / alpha);
    phase = 2 * pi * rand(M, numSymbols);
    H = magnitude .* exp(1j * phase);
end