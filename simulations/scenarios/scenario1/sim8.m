% Parameters
fc = 30e9;             % Carrier frequency (Hz, Ka-band)
c = 3e8;               % Speed of light (m/s)
Pt = 1;                % Transmit power (W)
Gt_dB = 35;            % Transmit antenna gain (dBi)
Gr_dB = 35;            % Receive antenna gain (dBi)
Gt = 10^(Gt_dB/10);
Gr = 10^(Gr_dB/10);
h = 600e3;             % Satellite altitude (m)
Re = 6371e3;           % Earth radius (m)
F = 10;                % Noise figure (linear)
k = 1.38e-23;          % Boltzmann constant
T = 290;               % Temperature (Kelvin)
BW = 100e6;            % Bandwidth (Hz)
N0 = k * T * F * BW;   % Noise power

% Doppler calculation
G = 6.67430e-11;       % Gravitational constant
M_E = 5.972e24;        % Mass of Earth
w_sat = sqrt(G*M_E/(Re + h)^3);
v_rel = w_sat * Re;    % Relative velocity at horizon
fd = fc * v_rel / c;   % Doppler shift (Hz)
fprintf('Estimated Doppler frequency: %.2f Hz\n', fd);

% Path loss calculation (free-space)
d = sqrt((Re + h)^2 - Re^2);
L_fs = (4 * pi * d * fc / c)^2;
Pr = Pt * Gt * Gr / L_fs;  % Received power (W)

% Simulation settings
numSymbols = 1e4;       % Number of symbols
K_dB = 10;              % Rician K-factor in dB
K = 10^(K_dB/10);
antenna_counts = [1 2 4 8]; % Number of receive antennas
avgBER_rayleigh = zeros(size(antenna_counts));
avgBER_rician = zeros(size(antenna_counts));

% Parameters for log-logistic fading
shape = 2;  % shape parameter alpha (adjust to change fading severity)
scale = 1;  % scale parameter beta

avgBER_loglogistic = zeros(size(antenna_counts));

for idx = 1:length(antenna_counts)
    M = antenna_counts(idx);
    
    % Rayleigh fading channel (M antennas)
    H_ray = (randn(M, numSymbols) + 1j*randn(M, numSymbols))/sqrt(2);
    H_MRC_ray = sum(abs(H_ray).^2,1);
    SNR_ray = (Pr .* H_MRC_ray) ./ (M * N0);
    ber_ray = 0.5 * erfc(sqrt(SNR_ray));
    avgBER_rayleigh(idx) = mean(ber_ray);
    
    % Rician fading channel (M antennas)
    LOS = sqrt(K/(K+1));
    NLOS = sqrt(1/(2*(K+1)));
    H_rician = LOS + NLOS*(randn(M,numSymbols) + 1j*randn(M,numSymbols));
    H_MRC_rician = sum(abs(H_rician).^2,1);
    SNR_rician = (Pr .* H_MRC_rician) ./ (M * N0);
    ber_rician = 0.5 * erfc(sqrt(SNR_rician));
    avgBER_rician(idx) = mean(ber_rician);
    
    % Log-logistic fading channel
    H_loglogistic = loglogistic_fading(M, numSymbols, shape, scale);
    H_MRC_loglogistic = sum(abs(H_loglogistic).^2,1);
    SNR_loglogistic = (Pr .* H_MRC_loglogistic) ./ (M * N0);
    ber_loglogistic = 0.5 * erfc(sqrt(SNR_loglogistic));
    avgBER_loglogistic(idx) = mean(ber_loglogistic);
end

% Plot results including log-logistic
figure;
semilogy(antenna_counts, avgBER_rayleigh, 'o-', 'LineWidth', 2);
hold on;
semilogy(antenna_counts, avgBER_rician, 's-', 'LineWidth', 2);
semilogy(antenna_counts, avgBER_loglogistic, 'd-', 'LineWidth', 2);
grid on;
xlabel('Number of Receive Antennas (M)');
ylabel('Average BER');
title('BER vs Number of Antennas (MRC, Ka-band, with Doppler)');
legend('Rayleigh', sprintf('Rician (K=%d dB)', K_dB), 'Log-logistic', 'Location', 'Best');
