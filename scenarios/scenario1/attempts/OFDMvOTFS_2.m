clc; clear; close all;

%% Parameters
N = 64;                % Number of subcarriers
M = 16;                % Number of OFDM symbols (time slots)
QAM_order = 4;         % QPSK (4-QAM)
SNR_dB = 0:5:30;       % SNR range for simulation

%% Channel parameters
fd_max = 200;          % Max Doppler shift (Hz) - can vary to simulate mobility
Ts = 1e-3;             % Symbol duration (1 ms)
fd_norm = fd_max * Ts; % Normalized Doppler frequency

path_delays = [0 1.5e-6 3.4e-6]; % seconds
path_gains_dB = [0 -8 -17];       % dB
path_gains = 10.^(path_gains_dB/20);

%% Generate QPSK symbols
M_qam = qammod(randi([0 QAM_order-1], N*M, 1), QAM_order, 'UnitAveragePower', true);

%% OFDM Modulation
ofdm_mod = @(x) ifft(reshape(x, N, M), N, 1) * sqrt(N);
ofdm_demod = @(y) fft(y, N, 1) / sqrt(N);

tx_ofdm = ofdm_mod(M_qam);

%% OTFS Modulation (basic implementation)
% OTFS modulates symbols in delay-Doppler domain using ISFFT
% ISFFT: Inverse Symplectic Finite Fourier Transform

% Define functions for OTFS
ISFFT = @(X) fft(ifft(X, N, 2), M, 1);
SFFT = @(Y) fft(ifft(Y, M, 1), N, 2);

% Arrange symbols in delay-Doppler grid
X_dd = reshape(M_qam, M, N).'; % N x M matrix

tx_otfs = ISFFT(X_dd);

%% Channel model: Multipath Rayleigh fading with Doppler
% Generate time-varying channel impulse response

% Sampling freq
Fs = N / Ts; % Depends on OFDM subcarrier spacing

t = (0:N*M-1) * Ts / N; % Time vector for N*M samples

% Generate channel impulse response h(t, tau)
h = zeros(length(path_delays), length(t));
for p = 1:length(path_delays)
    % Doppler shifts modeled as complex sinusoid with random phase
    fd_p = fd_max * randn * 0.5; % Random Doppler per path (simplified)
    h(p, :) = path_gains(p) * exp(1j*2*pi*fd_p*t);
end

% Construct channel impulse response in time domain (sum of delayed paths)
channel_response = zeros(1, length(t));
for p = 1:length(path_delays)
    delay_samples = round(path_delays(p)*Fs);
    h_p = [zeros(1,delay_samples), h(p,1:end-delay_samples)];
    channel_response = channel_response + h_p;
end

% Normalize channel power
channel_response = channel_response / norm(channel_response);

%% Transmit signal through channel
% OFDM
tx_ofdm_chan = conv(tx_ofdm(:).', channel_response, 'same');

% OTFS
tx_otfs_chan = conv(tx_otfs(:).', channel_response, 'same');

%% Add AWGN for different SNRs and compute BER
ber_ofdm = zeros(size(SNR_dB));
ber_otfs = zeros(size(SNR_dB));

for idx = 1:length(SNR_dB)
    noise_var = 10^(-SNR_dB(idx)/10);
    
    % Add noise
    rx_ofdm = tx_ofdm_chan + sqrt(noise_var/2)*(randn(size(tx_ofdm_chan)) + 1j*randn(size(tx_ofdm_chan)));
    rx_otfs = tx_otfs_chan + sqrt(noise_var/2)*(randn(size(tx_otfs_chan)) + 1j*randn(size(tx_otfs_chan)));
    
    % Receiver processing
    % OFDM demod
    rx_ofdm_mat = reshape(rx_ofdm, N, M);
    rx_ofdm_demod = ofdm_demod(rx_ofdm_mat);
    rx_ofdm_symbols = rx_ofdm_demod(:);
    
    % OTFS demod: apply SFFT to received symbols
    rx_otfs_mat = reshape(rx_otfs, N, M);
    rx_otfs_demod = SFFT(rx_otfs_mat);
    rx_otfs_symbols = rx_otfs_demod(:);
    
    % Demodulate QPSK
    rx_ofdm_data = qamdemod(rx_ofdm_symbols, QAM_order, 'UnitAveragePower', true);
    rx_otfs_data = qamdemod(rx_otfs_symbols, QAM_order, 'UnitAveragePower', true);
    
    % Compute BER
    ber_ofdm(idx) = mean(rx_ofdm_data ~= qamdemod(M_qam, QAM_order, 'UnitAveragePower', true));
    ber_otfs(idx) = mean(rx_otfs_data ~= qamdemod(M_qam, QAM_order, 'UnitAveragePower', true));
end

%% Plot results
figure;
semilogy(SNR_dB, ber_ofdm, '-o', 'LineWidth', 2);
hold on;
semilogy(SNR_dB, ber_otfs, '-s', 'LineWidth', 2);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
legend('OFDM', 'OTFS');
title(sprintf('BER vs SNR (Max Doppler = %d Hz)', fd_max));

