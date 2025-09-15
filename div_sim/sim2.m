clc; clear; close all;

% Parameters
N = 1e5;
SNR_dB = 0:2:20;
SNR = 10.^(SNR_dB/10);
M = 4; % QPSK
k = log2(M);
ber_nodiv = zeros(size(SNR));
ber_sc = zeros(size(SNR));
ber_mrc = zeros(size(SNR));

fd = 100;                 % Doppler frequency (Hz)
Ts = 1e-3;                % Symbol time
rho_doppler = besselj(0, 2*pi*fd*Ts);  % Doppler correlation

channel_corr = 0.7;       % Correlation coefficient between channels

% QPSK symbol mapping (Gray)
qpsk_sym = [1+1i, -1+1i, -1-1i, 1-1i] / sqrt(2);
bit2sym = @(b) qpsk_sym(1 + b(1)*2 + b(2));

% Generate QPSK symbols
bits = randi([0 1], 2, N);
symbols = zeros(1,N);
for i = 1:N
    symbols(i) = bit2sym(bits(:,i));
end

% Initial channels
h1 = (randn(1,N) + 1i*randn(1,N)) / sqrt(2);
h2 = channel_corr * h1 + sqrt(1 - channel_corr^2) * ((randn(1,N) + 1i*randn(1,N)) / sqrt(2));

for k_idx = 1:length(SNR)
    noise_var = 1/SNR(k_idx);

    % Time-varying fading (Jakes approx.)
    for t = 2:N
        h1(t) = rho_doppler * h1(t-1) + sqrt(1 - rho_doppler^2) * ((randn + 1i*randn)/sqrt(2));
        h2(t) = channel_corr * h1(t) + sqrt(1 - channel_corr^2) * ((randn + 1i*randn)/sqrt(2));
    end

    % Noise
    n1 = sqrt(noise_var/2) * (randn(1,N) + 1i*randn(1,N));
    n2 = sqrt(noise_var/2) * (randn(1,N) + 1i*randn(1,N));

    % Received signals
    r1 = h1 .* symbols + n1;
    r2 = h2 .* symbols + n2;

    % No diversity
    r_eq = r1 ./ h1;
    bits_rx = qpsk_demod(r_eq);
    ber_nodiv(k_idx) = sum(bits_rx(:) ~= bits(:)) / (2*N);

    % SC
    snr1 = abs(h1).^2;
    snr2 = abs(h2).^2;
    use_first = snr1 > snr2;
    r_sc = r1; r_sc(~use_first) = r2(~use_first);
    h_sc = h1; h_sc(~use_first) = h2(~use_first);
    r_eq_sc = r_sc ./ h_sc;
    bits_rx_sc = qpsk_demod(r_eq_sc);
    ber_sc(k_idx) = sum(bits_rx_sc(:) ~= bits(:)) / (2*N);

    % MRC
    r_mrc = conj(h1).*r1 + conj(h2).*r2;
    h_mrc = abs(h1).^2 + abs(h2).^2;
    r_eq_mrc = r_mrc ./ h_mrc;
    bits_rx_mrc = qpsk_demod(r_eq_mrc);
    ber_mrc(k_idx) = sum(bits_rx_mrc(:) ~= bits(:)) / (2*N);
end

% Plot
figure;
semilogy(SNR_dB, ber_nodiv, 'r-o', 'LineWidth', 2); hold on;
semilogy(SNR_dB, ber_sc, 'b-*', 'LineWidth', 2);
semilogy(SNR_dB, ber_mrc, 'g-s', 'LineWidth', 2);
xlabel('SNR (dB)'); ylabel('BER'); grid on;
legend('No Diversity', 'Selection Combining', 'MRC');
title('Satellite Diversity with QPSK, Doppler and Channel Correlation');

% ---- QPSK Demodulation Function ----
function bits_out = qpsk_demod(sym)
    re = real(sym); im = imag(sym);
    b1 = re > 0;
    b2 = im > 0;
    bits_out = [b1; b2];
end
