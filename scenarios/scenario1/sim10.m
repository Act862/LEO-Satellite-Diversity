% Constants
fc = 30e9;              % Carrier frequency (Hz, Ka-band)
c = 3e8;                % Speed of light (m/s)
Pt = 1;                 % Transmit power (W)
Gt_dB = 35;             % Transmit antenna gain (dBi)
Gr_dB = 35;             % Receive antenna gain (dBi)
Gt = 10^(Gt_dB/10);
Gr = 10^(Gr_dB/10);
Re = 6371e3;            % Earth radius (m)
h = 600e3;              % Satellite altitude (m)
F = 10;                 % Noise figure (linear)
k = 1.38e-23;           % Boltzmann constant
T = 290;                % Temperature (K)
BW = 100e6;             % Bandwidth (Hz)
N0 = k * T * F * BW;
K_dB = 10;              % Rician K-factor in dB
K = 10^(K_dB/10);

% Time vector simulating satellite pass (horizon -> overhead -> horizon)
time = linspace(0, 1, 200);         % Normalized pass time (0 to 1)
theta_deg = 90 * sin(pi * time);   % Elevation angle over time
theta_rad = deg2rad(theta_deg);

% Slant range function
slant_range = @(theta) sqrt((Re + h).^2 - (Re .* cos(theta)).^2) - Re .* sin(theta);
d = slant_range(theta_rad);

% Free-space path loss
Lfs = (4 * pi * d * fc / c).^2;
Pr = Pt * Gt * Gr ./ Lfs;

% Simulation parameters
M = 4;                  % Number of receive antennas
N = 10000;              % Fading samples
avgBER_rayleigh = zeros(size(time));
avgBER_rician = zeros(size(time));
SNR_rayleigh_dB = zeros(size(time));
SNR_rician_dB = zeros(size(time));

for idx = 1:length(time)
    pr = Pr(idx);

    % Rayleigh fading
    H_ray = (randn(M,N) + 1j*randn(M,N)) / sqrt(2);
    H_MRC_ray = sum(abs(H_ray).^2, 1);
    SNR_ray = (pr * H_MRC_ray) / (M * N0);
    SNR_rayleigh_dB(idx) = 10*log10(mean(SNR_ray));
    avgBER_rayleigh(idx) = mean(0.5 * erfc(sqrt(SNR_ray)));

    % Rician fading
    LOS = sqrt(K / (K + 1));
    NLOS = sqrt(1 / (2 * (K + 1)));
    H_rician = LOS + NLOS * (randn(M,N) + 1j*randn(M,N));
    H_MRC_rician = sum(abs(H_rician).^2, 1);
    SNR_rician = (pr * H_MRC_rician) / (M * N0);
    SNR_rician_dB(idx) = 10*log10(mean(SNR_rician));
    avgBER_rician(idx) = mean(0.5 * erfc(sqrt(SNR_rician)));
end

% Plotting SNR
figure;
plot(time, SNR_rayleigh_dB, 'b-', 'LineWidth', 2); hold on;
plot(time, SNR_rician_dB, 'r--', 'LineWidth', 2);
xlabel('Normalized Time (Satellite Pass)');
ylabel('Average SNR (dB)');
title('SNR vs Satellite Pass Time');
legend('Rayleigh', 'Rician (K=10dB)');
grid on;

% Plotting BER
figure;
semilogy(time, avgBER_rayleigh, 'b-', 'LineWidth', 2); hold on;
semilogy(time, avgBER_rician, 'r--', 'LineWidth', 2);
xlabel('Normalized Time (Satellite Pass)');
ylabel('Average BER');
title('BER vs Satellite Pass Time');
legend('Rayleigh', 'Rician (K=10dB)');
grid on;
