clc;
clear;
close all;

%% Parameters
fc = 2e9;                 % Carrier frequency (2 GHz)
c = 3e8;                  % Speed of light (m/s)
lambda = c / fc;          % Wavelength

v = 30;                   % Relative speed (m/s) for Doppler
theta_deg = -180:1:180;   % Angle range for Doppler spectrum

%% === Delay Spread Simulation ===
% Multipath parameters
num_paths = 5;
delays = [0 100 200 350 500]*1e-9; % delays in seconds (ns scale)
powers_dB = [0 -3 -6 -9 -12];      % path power levels in dB
powers = 10.^(powers_dB/10);      % convert to linear scale

% Compute mean delay and RMS delay spread
mean_delay = sum(powers .* delays) / sum(powers);
mean_sq_delay = sum(powers .* delays.^2) / sum(powers);
tau_rms = sqrt(mean_sq_delay - mean_delay^2);

% Plot Power Delay Profile (PDP)
figure;
stem(delays*1e9, powers, 'filled');
title(['Power Delay Profile — RMS Delay Spread: ' num2str(tau_rms*1e9, '%.2f') ' ns']);
xlabel('Delay (ns)');
ylabel('Normalized Power');
grid on;

%% === Doppler Spread Simulation ===
% Doppler shift as function of angle
theta = deg2rad(theta_deg);
f_D = (v / c) * fc * cos(theta);

% Doppler Power Spectrum (Jakes' Model)
S_f = 1 ./ sqrt(1 - (f_D / max(abs(f_D))).^2);
S_f(~isfinite(S_f)) = 0; % clean up infinities

% Plot Doppler Spectrum
figure;
plot(f_D, S_f, 'b', 'LineWidth', 2);
title(['Doppler Power Spectrum — Max Doppler Shift: ' num2str((v/c)*fc, '%.2f') ' Hz']);
xlabel('Doppler Shift (Hz)');
ylabel('Spectral Density (a.u.)');
xlim([-max(abs(f_D)) max(abs(f_D))]);
ylim([0 max(S_f)*1.1]);
grid on;

