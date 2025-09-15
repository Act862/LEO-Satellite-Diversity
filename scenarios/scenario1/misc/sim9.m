% Constants
fc = 30e9;               % Carrier frequency (Hz)
c = 3e8;                 % Speed of light (m/s)
Pt = 1;                  % Transmit power (W)
Gt_dB = 35;              % Transmit gain (dBi)
Gr_dB = 35;              % Receive gain (dBi)
Gt = 10^(Gt_dB/10);
Gr = 10^(Gr_dB/10);
Re = 6371e3;             % Earth radius (m)
h = 600e3;               % Satellite altitude (m)
F = 10;                  % Noise figure (linear)
k = 1.38e-23;            % Boltzmann constant
T = 290;                 % System temperature (K)
BW = 100e6;              % Bandwidth (Hz)
N0 = k * T * F * BW;     % Noise power
K_dB = 10;               % Rician K-factor in dB
K = 10^(K_dB/10);

% Elevation angles
theta_deg = linspace(5, 90, 100);
theta_rad = deg2rad(theta_deg);

% Slant range function
slant_range = @(theta) sqrt((Re + h).^2 - (Re .* cos(theta)).^2) - Re .* sin(theta);
d = slant_range(theta_rad);

% Free-space path loss
L_fs = (4 * pi * d * fc / c).^2;
Pr = Pt * Gt * Gr ./ L_fs;

% Simulation parameters
M = 4;                  % Number of antennas
N = 1e4;                % Symbols per elevation angle
avgBER_rayleigh = zeros(size(theta_deg));
avgBER_rician = zeros(size(theta_deg));

for i = 1:length(theta_deg)
    pr = Pr(i);
    
    % Rayleigh fading
    H_ray = (randn(M, N) + 1j*randn(M, N)) / sqrt(2);
    H_MRC_ray = sum(abs(H_ray).^2, 1);
    SNR_ray = (pr * H_MRC_ray) / (M * N0);
    ber_ray = 0.5 * erfc(sqrt(SNR_ray));
    avgBER_rayleigh(i) = mean(ber_ray);

    % Rician fading
    LOS = sqrt(K / (K + 1));
    NLOS = sqrt(1 / (2 * (K + 1)));
    H_rician = LOS + NLOS * (randn(M, N) + 1j * randn(M, N));
    H_MRC_rician = sum(abs(H_rician).^2, 1);
    SNR_rician = (pr * H_MRC_rician) / (M * N0);
    ber_rician = 0.5 * erfc(sqrt(SNR_rician));
    avgBER_rician(i) = mean(ber_rician);
end

% Plotting
figure;
semilogy(theta_deg, avgBER_rayleigh, 'b-', 'LineWidth', 2); hold on;
semilogy(theta_deg, avgBER_rician, 'r--', 'LineWidth', 2);
xlabel('Elevation Angle (degrees)');
ylabel('Average BER');
title('BER vs Elevation Angle (Ka-band, LEO, MRC with 4 antennas)');
legend('Rayleigh', sprintf('Rician (K = %d dB)', K_dB), 'Location', 'SouthWest');
grid on;
