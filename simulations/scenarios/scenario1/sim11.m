% Constants
fc = 30e9;               % Carrier frequency (Hz, Ka-band)
c = 3e8;                 % Speed of light (m/s)
Pt = 1;                  % Transmit power (W)
Gt_dB = 35;              % Transmit antenna gain (dBi)
Gr_dB = 35;              % Receive antenna gain (dBi)
Gt = 10^(Gt_dB/10);
Gr = 10^(Gr_dB/10);
Re = 6371e3;             % Earth radius (m)
h_sat = 600e3;           % Satellite altitude (m)
F = 10;                  % Noise figure (linear)
k = 1.38e-23;            % Boltzmann constant
T = 290;                 % Temperature (K)
BW = 100e6;              % Bandwidth (Hz)
N0 = k * T * F * BW;
K_dB = 10;               % Rician K-factor in dB
K = 10^(K_dB/10);

% Altitudes of VSATs
h_VSAT_low = 0;          % Sea level
h_VSAT_high = 3000;      % 3000 meters above sea level

% Satellite pass (normalized time)
time = linspace(0, 1, 200);
theta_deg = 90 * sin(pi * time);   % Elevation angle (0 to 90 to 0)
theta_rad = deg2rad(theta_deg);

% Function to compute slant range for a given VSAT altitude
slant_range = @(theta, h_vsat) sqrt((Re + h_sat).^2 - (Re + h_vsat).^2 .* (cos(theta)).^2) ...
                                - (Re + h_vsat) .* sin(theta);

% Compute slant ranges
d_low = slant_range(theta_rad, h_VSAT_low);
d_high = slant_range(theta_rad, h_VSAT_high);

% Compute path loss and received power
Lfs_low = (4 * pi * d_low * fc / c).^2;
Lfs_high = (4 * pi * d_high * fc / c).^2;
Pr_low = Pt * Gt * Gr ./ Lfs_low;
Pr_high = Pt * Gt * Gr ./ Lfs_high;

% Fading model
M = 4; N = 10000;
avgBER_low = zeros(size(time));
avgBER_high = zeros(size(time));
SNR_low_dB = zeros(size(time));
SNR_high_dB = zeros(size(time));

for idx = 1:length(time)
    % Low-altitude VSAT
    pr1 = Pr_low(idx);
    LOS = sqrt(K / (K + 1));
    NLOS = sqrt(1 / (2 * (K + 1)));
    H1 = LOS + NLOS * (randn(M,N) + 1j*randn(M,N));
    Hsum1 = sum(abs(H1).^2, 1);
    SNR1 = pr1 * Hsum1 / (M * N0);
    avgBER_low(idx) = mean(0.5 * erfc(sqrt(SNR1)));
    SNR_low_dB(idx) = 10 * log10(mean(SNR1));

    % High-altitude VSAT
    pr2 = Pr_high(idx);
    H2 = LOS + NLOS * (randn(M,N) + 1j*randn(M,N));
    Hsum2 = sum(abs(H2).^2, 1);
    SNR2 = pr2 * Hsum2 / (M * N0);
    avgBER_high(idx) = mean(0.5 * erfc(sqrt(SNR2)));
    SNR_high_dB(idx) = 10 * log10(mean(SNR2));
end

% Plot SNR Comparison
figure;
plot(time, SNR_low_dB, 'b-', 'LineWidth', 2); hold on;
plot(time, SNR_high_dB, 'r--', 'LineWidth', 2);
xlabel('Normalized Time (Satellite Pass)');
ylabel('Average SNR (dB)');
title('SNR vs Time: Sea Level vs Mountain VSAT');
legend('VSAT at Sea Level', 'VSAT at 3000 m Altitude');
grid on;

% Plot BER Comparison
figure;
semilogy(time, avgBER_low, 'b-', 'LineWidth', 2); hold on;
semilogy(time, avgBER_high, 'r--', 'LineWidth', 2);
xlabel('Normalized Time (Satellite Pass)');
ylabel('Average BER');
title('BER vs Time: Sea Level vs Mountain VSAT');
legend('VSAT at Sea Level', 'VSAT at 3000 m Altitude');
grid on;
