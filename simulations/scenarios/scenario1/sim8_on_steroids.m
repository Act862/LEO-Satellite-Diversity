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
F_dB = 10;             % Noise figure (dB)
F = 10^(F_dB/10);      % Noise figure (linear)
k = 1.38e-23;          % Boltzmann constant
T = 290;               % Temperature (Kelvin)
BW = 100e6;            % Bandwidth (Hz)

% Elevation angle and slant range
el_deg = 90;           % Elevation angle in degrees
el_rad = deg2rad(el_deg);
d = sqrt((Re + h)^2 - (Re * cos(el_rad))^2) - Re * sin(el_rad); % Slant range

% Doppler calculation
G = 6.67430e-11;       % Gravitational constant
M_E = 5.972e24;        % Mass of Earth
w_sat = sqrt(G*M_E/(Re + h)^3);
v_rel = w_sat * Re;    % Relative velocity at horizon
fd = fc * v_rel / c;   % Doppler shift (Hz)
fprintf('Estimated Doppler frequency: %.2f Hz\n', fd);

% Path loss calculation (free-space)
L_fs = (4 * pi * d * fc / c)^2;

% Total received power with all losses
Pr = Pt * Gt * Gr / L_fs;  % Received power (W)

% Noise power
N0 = k * T * F * BW;   % Noise power

% Simulation settings
numSymbols = 1e4;       % Number of symbols
K_dB = 10;              % Rician K-factor in dB
K = 10^(K_dB/10);
antenna_counts = [1 2 4 8]; % Number of receive antennas
avgBER_rayleigh = zeros(size(antenna_counts));
avgBER_rician = zeros(size(antenna_counts));
avgBER_loglogistic = zeros(size(antenna_counts));
avgOutage = zeros(size(antenna_counts));
outage_threshold_dB = 5; % SNR outage threshold in dB
outage_threshold_linear = 10^(outage_threshold_dB / 10);

% Parameters for log-logistic fading
shape = 2;  % shape parameter alpha
scale = 1;  % scale parameter beta

for idx = 1:length(antenna_counts)
    M = antenna_counts(idx);

    % Rayleigh fading channel with Doppler
    H_ray = doppler_fading(numSymbols, fd, M);
    H_MRC_ray = sum(abs(H_ray).^2, 1);
    SNR_ray = (Pr .* H_MRC_ray) ./ (M * N0);
    ber_ray = 0.5 * erfc(sqrt(SNR_ray));
    avgBER_rayleigh(idx) = mean(ber_ray);

    % Outage calculation
    avgOutage(idx) = mean(10*log10(SNR_ray) < outage_threshold_dB);

    % Rician fading channel with Doppler
    H_rician = doppler_fading(numSymbols, fd, M, K);
    H_MRC_rician = sum(abs(H_rician).^2, 1);
    SNR_rician = (Pr .* H_MRC_rician) ./ (M * N0);
    ber_rician = 0.5 * erfc(sqrt(SNR_rician));
    avgBER_rician(idx) = mean(ber_rician);

    % Log-logistic fading channel
    H_loglogistic = loglogistic_fading(M, numSymbols, shape, scale);
    H_MRC_loglogistic = sum(abs(H_loglogistic).^2, 1);
    SNR_loglogistic = (Pr .* H_MRC_loglogistic) ./ (M * N0);
    ber_loglogistic = 0.5 * erfc(sqrt(SNR_loglogistic));
    avgBER_loglogistic(idx) = mean(ber_loglogistic);
end

% Plot BER results
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

% Plot Outage Probability
figure;
plot(antenna_counts, avgOutage, 'p-', 'LineWidth', 2);
grid on;
xlabel('Number of Receive Antennas (M)');
ylabel('Outage Probability');
title(sprintf('Outage Probability vs Number of Antennas (Threshold = %d dB)', outage_threshold_dB));